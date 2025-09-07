# worker.py
import os, time, math, json
from datetime import datetime, timezone
from pymongo import MongoClient, ASCENDING, DESCENDING, errors

try:
    import openseespy.opensees as ops
    HAS_OPS = True
except Exception:
    HAS_OPS = False

# ===== Config =====
MONGO_URI   = os.getenv("MONGODB_URI", "")
MONGO_DB    = os.getenv("MONGODB_DB", "biostrucx")
CLIENTS     = [c.strip() for c in os.getenv("FEM_CLIENTS", "jeimie").split(",") if c.strip()]
INTERVAL_S  = int(os.getenv("FEM_INTERVAL_SEC", "120"))

# Geometría / propiedades (m, Pa, N/m)
L = float(os.getenv("FEM_L", "25"))
B = float(os.getenv("FEM_B", "0.25"))      # ancho (eje y)
H = float(os.getenv("FEM_H", "1.0"))       # “espesor” visual (altura)
E = float(os.getenv("FEM_E", str(30e9)))
rho = 2500.0
q_kNpm = float(os.getenv("FEM_W", "15"))
q = q_kNpm * 1e3

# Discretización
N     = int(os.getenv("FEM_N", "80"))      # a lo largo (x) -> N segmentos, N+1 nodos
Wdiv  = int(os.getenv("FEM_WDIV", "8"))    # a lo ancho (y)
Hdiv  = int(os.getenv("FEM_HDIV", "4"))    # a lo alto (espesor “vertical” alrededor de la deformada)
SENSOR_X_FRACTION = float(os.getenv("FEM_SENSOR_X", "0.65"))

# Escalas visuales
DEF_SCALE      = float(os.getenv("FEM_DEF_SCALE", "30.0"))  # amplifica la deformada
THICK_SCALE    = float(os.getenv("FEM_THICK_SCALE", "1.0")) # factor visual del espesor H (1=real)

def mongo():
    cli = MongoClient(MONGO_URI, serverSelectionTimeoutMS=8000)
    return cli[MONGO_DB]

def ensure_indexes(db):
    for name in ("simulation_result","simulation_ts"):
        try:
            db[name].create_index([("clientid", ASCENDING), ("ts", DESCENDING)])
        except errors.OperationFailure:
            pass

# ---- Viga simplemente apoyada con carga distribuida q ----
def deflection_analytic(x):
    I = (B * (H**3)) / 12.0
    return (q * x * (L**3 - 2*L*(x**2) + x**3)) / (24.0 * E * I)

def deflection_opensees():
    if not HAS_OPS:
        return None
    try:
        ndm, ndf = 2, 3
        ops.wipe()
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)

        # nodos a lo largo de la viga
        for i in range(N + 1):
            x = L * i / N
            ops.node(i + 1, x, 0.0)

        # apoyos simples (u=v=0, rot libre)
        ops.fix(1,     1, 1, 0)
        ops.fix(N + 1, 1, 1, 0)

        # sección equivalente y propiedades
        A = B * H
        I = (B * (H**3)) / 12.0

        # *** ESTA LÍNEA FALTABA ***
        # Transformation para elementos en 2D (Lineal o PDelta)
        ops.geomTransf('Linear', 1)

        # elementos viga elástica (usa el transfTag = 1)
        for i in range(1, N + 1):
            ops.element('elasticBeamColumn', i, i, i + 1, A, E, I, 1)

        # carga distribuida (aprox. con puntuales nodales)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        p = -q * (L / N)  # N por nodo (hacia abajo)
        for i in range(2, N):  # evita cargar apoyos
            ops.load(i, 0.0, p, 0.0)

        # análisis estático lineal
        ops.system('BandGeneral')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        ops.analyze(1)

        # deflexión vertical uy en cada nodo
        w_def = [ops.nodeDisp(i + 1, 2) for i in range(N + 1)]
        return w_def

    except Exception:
        # si algo falla, deja que el caller use el analítico
        return None


def build_viz():
    # 1) deflexión eje neutro
    w_def = deflection_opensees()
    if w_def is None:
        w_def = [deflection_analytic(L * i / N) for i in range(N + 1)]

    # 2) malla volumétrica tipo prisma rectangular:
    #    ejes: x (largo), y (ancho B), z (vertical = deformada + espesor)
    Ny = Wdiv + 1
    Nz = Hdiv + 1
    stride_x = Ny * Nz

    def vidx(i, j, k):
        return i * stride_x + j * Nz + k

    vertices = []
    u_mag = []

    # pre-cálculos
    halfB = B * 0.5
    halfH = (H * THICK_SCALE) * 0.5

    for i in range(N + 1):
        x = L * i / N
        w = w_def[i] * DEF_SCALE              # deformada visual
        mm = abs(w_def[i]) * 1000.0           # mm reales (para color)

        for j in range(Ny):                   # ancho
            y = -halfB + (B * j / Wdiv)
            for k in range(Nz):               # espesor “vertical” alrededor de w
                h = -halfH + (2 * halfH * k / Hdiv)
                z = w + h                     # top/bottom alrededor de la deformada
                vertices.extend([x, y, z])
                u_mag.append(mm)

    indices = []

    def add_quad(a, b, c, d):
        # dos triángulos: a-b-c y a-c-d
        indices.extend([a, b, c, a, c, d])

    # 3) carcasas:
    # top (k=Hdiv) y bottom (k=0)
    for i in range(N):
        for j in range(Wdiv):
            # top
            a = vidx(i,     j,     Hdiv)
            b = vidx(i + 1, j,     Hdiv)
            c = vidx(i + 1, j + 1, Hdiv)
            d = vidx(i,     j + 1, Hdiv)
            add_quad(a, b, c, d)
            # bottom
            a = vidx(i,     j,     0)
            b = vidx(i + 1, j,     0)
            c = vidx(i + 1, j + 1, 0)
            d = vidx(i,     j + 1, 0)
            add_quad(d, c, b, a)  # invertido para normal opuesta

    # lados y = -B/2 y y = +B/2
    for i in range(N):
        for k in range(Hdiv):
            # lado -Y (j=0)
            a = vidx(i,     0, k)
            b = vidx(i + 1, 0, k)
            c = vidx(i + 1, 0, k + 1)
            d = vidx(i,     0, k + 1)
            add_quad(d, c, b, a)
            # lado +Y (j=Wdiv)
            a = vidx(i,     Wdiv, k)
            b = vidx(i + 1, Wdiv, k)
            c = vidx(i + 1, Wdiv, k + 1)
            d = vidx(i,     Wdiv, k + 1)
            add_quad(a, b, c, d)

    # tapas (opcional): i=0 y i=N para cerrar extremos
    for j in range(Wdiv):
        for k in range(Hdiv):
            # i=0
            a = vidx(0, j,     k)
            b = vidx(0, j + 1, k)
            c = vidx(0, j + 1, k + 1)
            d = vidx(0, j,     k + 1)
            add_quad(a, b, c, d)
            # i=N
            a = vidx(N, j,     k)
            b = vidx(N, j + 1, k)
            c = vidx(N, j + 1, k + 1)
            d = vidx(N, j,     k + 1)
            add_quad(d, c, b, a)

    # 4) marcador (centro del espesor)
    xs = SENSOR_X_FRACTION * L
    k_near = min(range(N + 1), key=lambda ii: abs(L * ii / N - xs))
    z_marker = w_def[k_near] * DEF_SCALE
    marker = [xs, 0.0, z_marker]

    max_mm = max(abs(v) for v in w_def) * 1000.0

    return vertices, indices, u_mag, marker, max_mm

def run_once(db, clientid):
    now = datetime.now(timezone.utc)
    vertices, indices, u_mag, marker, max_mm = build_viz()

    db["simulation_result"].insert_one({
        "clientid": clientid,
        "status": "done",
        "ts": now,
        "model": { "type": "beam_demo", "L_m": L, "B_m": B, "H_m": H, "E_Pa": E, "q_Npm": q },
        "params": {},
        "viz": { "vertices": vertices, "indices": indices, "u_mag": u_mag, "marker": marker }
    })

    db["simulation_ts"].insert_one({ "ts": now, "clientid": clientid, "fem_mm": float(max_mm) })

def main():
    db = mongo()
    ensure_indexes(db)
    print(f"[worker] starting… clients={CLIENTS} every {INTERVAL_S}s")
    while True:
        for cid in CLIENTS:
            try:
                run_once(db, cid)
            except Exception as e:
                print("[worker] error:", e)
        time.sleep(INTERVAL_S)

if __name__ == "__main__":
    main()

