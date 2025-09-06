import os, time, math
from datetime import datetime, timezone
from pymongo import MongoClient, ASCENDING

# ---------- Env ----------
MONGODB_URI = os.environ.get("MONGODB_URI")
MONGODB_DB  = os.environ.get("MONGODB_DB", "biostrucx")
FEM_CLIENTS = [c.strip() for c in os.environ.get("FEM_CLIENTS", "jeimie").split(",") if c.strip()]
INTERVAL    = int(os.environ.get("FEM_INTERVAL_SEC", "120"))

# Parámetros del demo (puedes moverlos a env si quieres)
L  = float(os.environ.get("FEM_L", 25.0))     # m
B  = float(os.environ.get("FEM_B", 0.25))     # m (ancho)
H  = float(os.environ.get("FEM_H", 0.25))     # m (alto)
E  = float(os.environ.get("FEM_E", 30e9))     # Pa
rho= float(os.environ.get("FEM_RHO", 2500.0)) # kg/m3 (no imprescindible aquí)
qkN= float(os.environ.get("FEM_W_KNM", 15.0)) # kN/m (carga distribuida)
N  = int(os.environ.get("FEM_N", 40))         # nº de divisiones (N+1 nodos)

# Escala de la deformada para que "se vea"
PLOT_SCALE = float(os.environ.get("FEM_PLOT_SCALE", 100.0))  # multiplica los desplazamientos (m) para dibujar

# Intentamos cargar OpenSeesPy
USE_OPS = True
try:
    from openseespy import opensees as ops
except Exception as e:
    print("[worker] OpenSeesPy no disponible, usaré fórmula analítica. Motivo:", e)
    USE_OPS = False

# ---------- Helpers de malla (ribbon triangulado) ----------
def ribbon_from_line(xs, ys, halfW=0.12):
    """
    Construye una malla tipo ribbon (cinta) extruyendo +/- halfW en z.
    xs, ys: listas del eje de la viga y su deformada (en m).
    Devuelve (vertices_flat, indices, u_mag_mm)
    """
    V = []
    U = []
    for x, y in zip(xs, ys):
        # dos vértices por estación (± en z)
        V.extend([x, y, -halfW])
        V.extend([x, y,  halfW])
        # magnitud para color (mm)
        umm = abs(y) * 1000.0
        U.extend([umm, umm])

    # Triangulación: (i0,i2,i1) y (i1,i2,i3) para cada panel
    idx = []
    for i in range(0, (len(xs) - 1) * 2, 2):
        i0, i1, i2, i3 = i, i + 1, i + 2, i + 3
        idx.extend([i0, i2, i1,  i1, i2, i3])

    return V, idx, U

# ---------- Solución analítica (fallback) ----------
def analytical_beam(L, q, E, I, N):
    """
    Viga simplemente apoyada, carga uniforme q (N/m), deflexión vertical w(x) (m).
    w(x) = q x (L^3 - 2 L x^2 + x^3) / (24 E I)
    """
    xs = [L * i / N for i in range(N + 1)]
    ws = []
    for x in xs:
        w = q * x * (L**3 - 2*L*(x**2) + x**3) / (24.0 * E * I)  # m
        ws.append(w)
    return xs, ws

# ---------- Solución con OpenSeesPy 2D y malla ----------
def solve_opensees_2d(L, B, H, E, q, N, plot_scale):
    """
    Modelo 2D (ndm=2, ndf=3) con elasticBeamColumn y carga uniforme.
    Devuelve xs (m), w_deformed (m, amplificada ya para dibujar) y w_real (m, real sin escala)
    """
    A = B * H
    I = (B * H**3) / 12.0  # flexión en 2D

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)  # 2D: Ux, Uy, Rz

    # Nodos sobre eje x, y=0
    node_tags = []
    for i in range(N + 1):
        x = L * i / N
        tag = i + 1
        ops.node(tag, x, 0.0)
        node_tags.append(tag)

    # Apoyos: simple (pin-roller)
    # i=0 (izquierda): pin -> Ux=Uy=0
    ops.fix(node_tags[0], 1, 1, 0)
    # i=N (derecha): roller -> Uy=0
    ops.fix(node_tags[-1], 0, 1, 0)

    # Geometría y elemento
    ops.geomTransf('Linear', 1)
    for i in range(N):
        ops.element('elasticBeamColumn', i + 1, node_tags[i], node_tags[i+1], A, E, I, 1)

    # Carga uniforme vertical (negativa hacia abajo en 2D: eje y)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    # en 2D: -beamUniform wy  (OpenSeesPy usa wy>0 hacia -Y; ponemos q positivo y signo "-" aquí)
    for i in range(1, N + 1):
        ops.eleLoad('-ele', i, '-type', '-beamUniform', q)  # q positivo = hacia -Y

    # Análisis estático
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("analyze() no convergió")

    # Coords y desplazamientos
    xs  = []
    w_r = []   # deformada real (m)
    for tag in node_tags:
        x, y = ops.nodeCoord(tag)
        ux, uy, rz = ops.nodeDisp(tag)
        xs.append(x)
        w_r.append(uy)   # vertical

    # para dibujar, amplificamos
    w_d = [plot_scale * wy for wy in w_r]

    return xs, w_d, w_r

# ---------- Core: resuelve y publica ----------
def run_once(db, clientid):
    col_res = db["simulation_result"]
    col_ts  = db["simulation_ts"]

    # Unidades
    q = qkN * 1000.0   # kN/m -> N/m
    A = B * H
    I = (B * H**3) / 12.0

    try:
        if USE_OPS:
            xs, w_draw, w_real = solve_opensees_2d(L, B, H, E, q, N, PLOT_SCALE)
        else:
            xs, w_real = analytical_beam(L, q, E, I, N)
            w_draw = [PLOT_SCALE * wy for wy in w_real]

        # Viga "ribbon" triangulada para tu viewer
        vertices, indices, u_mag = ribbon_from_line(xs, w_draw, halfW=B*0.5)

        # Deflexión en el centro (real, m → mm)
        mid_idx = N // 2
        fem_mm = abs(w_real[mid_idx]) * 1000.0

        # Marcador: por ejemplo a 0.8*L
        mark_x  = 0.8 * L
        # interpola deformada para ese x
        j = min(N-1, int((mark_x / L) * N))
        xa, xb = xs[j], xs[j+1]
        ya, yb = w_draw[j], w_draw[j+1]
        t = 0.0 if xb == xa else (mark_x - xa) / (xb - xa)
        mark_y = ya*(1-t) + yb*t
        marker  = [mark_x, mark_y, 0.0]

        doc = {
            "clientid": clientid,
            "ts": datetime.now(timezone.utc),
            "status": "done",
            "model": {
                "type": "beam_2d_linear",
                "params": {"L": L, "B": B, "H": H, "E": E, "q_Npm": q, "N": N}
            },
            "viz": {
                "vertices": vertices,   # [x,y,z,...]
                "indices":  indices,    # [i0,i1,i2,...]
                "u_mag":    u_mag,      # [mm,...] (para color)
                "marker":   marker
            }
        }

        # Upsert del resultado “latest”
        col_res.update_one(
            {"clientid": clientid},
            {"$set": doc},
            upsert=True
        )

        # Serie temporal (para el gráfico 1)
        col_ts.insert_one({
            "clientid": clientid,
            "ts": datetime.now(timezone.utc),
            "fem_mm": fem_mm
        })

        print(f"[worker] {clientid} ok | midspan = {fem_mm:.3f} mm")

    except Exception as e:
        print(f"[worker] {clientid} ERROR:", e)

def main():
    client = MongoClient(MONGODB_URI)
    db = client[MONGODB_DB]
    db["simulation_ts"].create_index([("clientid", ASCENDING), ("ts", ASCENDING)])

    print(f"[worker] starting… clients={FEM_CLIENTS} every {INTERVAL} sec")
    while True:
        for cid in FEM_CLIENTS:
            run_once(db, cid)
        time.sleep(INTERVAL)

if __name__ == "__main__":
    main()

