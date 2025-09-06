# worker.py
import os, time, math, json
from datetime import datetime, timezone
from pymongo import MongoClient, ASCENDING, DESCENDING, errors

try:
    # Intentar OpenSeesPy; si no está, caemos al analítico
    import openseespy.opensees as ops
    HAS_OPS = True
except Exception:
    HAS_OPS = False

# ===== Parámetros por defecto (puedes sobreescribir con ENV) =====
MONGO_URI   = os.getenv("MONGODB_URI", "")
MONGO_DB    = os.getenv("MONGODB_DB", "biostrucx")
CLIENTS     = [c.strip() for c in os.getenv("FEM_CLIENTS", "jeimie").split(",") if c.strip()]
INTERVAL_S  = int(os.getenv("FEM_INTERVAL_SEC", "120"))

# Geometría/propiedades del demo (m)
L = float(os.getenv("FEM_L", "25"))      # largo
B = float(os.getenv("FEM_B", "0.25"))    # ancho “visual” del ribbon (no sección)
H = float(os.getenv("FEM_H", "1.0"))     # alto (para I)
E = float(os.getenv("FEM_E", str(30e9))) # Pa
rho = 2500.0
q_kNpm = float(os.getenv("FEM_W", "15")) # carga distribuida (kN/m)
q = q_kNpm * 1e3                         # N/m

# Malla de visualización
N = int(os.getenv("FEM_N", "80"))        # divisiones a lo largo (>= 2)
Wdiv = int(os.getenv("FEM_WDIV", "8"))   # divisiones a lo ancho (>= 2)
SENSOR_X_FRACTION = float(os.getenv("FEM_SENSOR_X", "0.65")) # x del marcador

# Escala de deformada para que “se vea”
DEF_SCALE = float(os.getenv("FEM_DEF_SCALE", "30.0"))

def mongo():
    cli = MongoClient(MONGO_URI, serverSelectionTimeoutMS=8000)
    db = cli[MONGO_DB]
    return db

def ensure_indexes(db):
    # Evita chocar con índices ya existentes
    try:
        db["simulation_result"].create_index([("clientid", ASCENDING), ("ts", DESCENDING)])
    except errors.OperationFailure:
        pass
    try:
        db["simulation_ts"].create_index([("clientid", ASCENDING), ("ts", DESCENDING)])
    except errors.OperationFailure:
        pass

# ====== Modelo simple: viga simplemente apoyada con q ======
def deflection_analytic(x):
    # w(x) en metros (teoría vigas Euler-Bernoulli)
    I = (B * (H**3)) / 12.0
    # fórmula: q*x*(L^3 - 2*L*x^2 + x^3) / (24*E*I)
    return (q * x * (L**3 - 2*L*(x**2) + x**3)) / (24.0 * E * I)

def deflection_opensees():
    # Discretización lineal para sacar forma modal estática
    if not HAS_OPS:
        return None
    try:
        ndm, ndf = 2, 3
        ops.wipe()
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)

        # Nnod = N+1
        for i in range(N + 1):
            x = L * i / N
            ops.node(i + 1, x, 0.0)

        # apoyos simples en 1 y N+1
        ops.fix(1, 1, 1, 0)            # empotramiento parcial (u=v=0, libre rot)
        ops.fix(N + 1, 1, 1, 0)

        # sección equivalente (elasticBeamColumn 2D)
        A = B * H
        I = (B * (H**3)) / 12.0
        ops.uniaxialMaterial('Elastic', 1, E)
        ops.section('Elastic', 1, E, A, I)

        # Elementos
        for i in range(1, N + 1):
            ops.element('elasticBeamColumn', i, i, i + 1, A, E, I, 1)

        # Carga distribuida (equivalente nodal)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        w = -q  # hacia abajo
        # repartir como puntuales: cada tramo L/N con carga w*(L/N)
        p = w * (L / N)
        for i in range(2, N):  # no cargar apoyos
            ops.load(i, 0.0, p, 0.0)

        ops.system('BandGeneral')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        ops.analyze(1)

        w_def = []
        for i in range(N + 1):
            ntag = i + 1
            ux = ops.nodeDisp(ntag, 1)
            uy = ops.nodeDisp(ntag, 2)
            # usamos uy (vertical) como “deflexión” (metros)
            w_def.append(uy)
        return w_def
    except Exception:
        return None

def build_viz():
    # 1) deflexión a lo largo
    w_def = deflection_opensees()
    if w_def is None:
        # analítico
        w_def = [deflection_analytic(L * i / N) for i in range(N + 1)]

    # 2) construir malla ribbon N×Wdiv
    vertices = []
    u_mag = []
    for i in range(N + 1):
        x = L * i / N
        w = w_def[i] * DEF_SCALE            # amplificado para verlo
        for j in range(Wdiv + 1):
            t = j / Wdiv                    # 0..1
            y = (t - 0.5) * B               # centrado en 0
            z = w                           # deformada en z
            vertices.extend([x, y, z])
            u_mag.append(abs(w_def[i] * 1000.0))  # magnitud real en mm (para color)

    # Triangulación
    indices = []
    stride = Wdiv + 1
    for i in range(N):
        for j in range(Wdiv):
            a = i * stride + j
            b = a + 1
            c = (i + 1) * stride + j
            d = c + 1
            # dos triángulos por quad
            indices.extend([a, c, b])
            indices.extend([b, c, d])

    # Marcador (posición del sensor)
    xs = SENSOR_X_FRACTION * L
    # buscamos z en esa x:
    k = min(range(N + 1), key=lambda ii: abs(L * ii / N - xs))
    z_marker = w_def[k] * DEF_SCALE
    marker = [xs, 0.0, z_marker]

    return vertices, indices, u_mag, marker, max(abs(v) for v in w_def) * 1000.0 # max mm

def run_once(db, clientid):
    now = datetime.now(timezone.utc)

    vertices, indices, u_mag, marker, max_mm = build_viz()

    # documento principal (para el visor)
    doc = {
        "clientid": clientid,
        "status": "done",
        "ts": now,
        "model": {
            "type": "beam_demo",
            "L_m": L, "B_m": B, "H_m": H, "E_Pa": E, "q_Npm": q
        },
        "params": {},
        "viz": {
            "vertices": vertices,   # float[]
            "indices": indices,     # int[] (triángulos)
            "u_mag": u_mag,         # float[] (mm por vértice)
            "marker": marker        # [x,y,z]
        }
    }
    db["simulation_result"].insert_one(doc)

    # punto de serie temporal (máx deflexión en mm)
    db["simulation_ts"].insert_one({
        "ts": now, "clientid": clientid, "fem_mm": float(max_mm)
    })

def main():
    db = mongo()
    ensure_indexes(db)

    print(f"[worker] starting… clients={CLIENTS} every {INTERVAL_S} sec")
    while True:
        for cid in CLIENTS:
            try:
                run_once(db, cid)
            except Exception as e:
                print("[worker] error:", e)
        time.sleep(INTERVAL_S)

if __name__ == "__main__":
    main()
