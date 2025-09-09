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

# =============================
# Geometría / propiedades (m, Pa, N/m)
# 👉 Ahora valores fijos (ignora ENV)
# =============================
L = 1.0       # longitud de la viga (m)
B = 0.30       # ancho (m)  ← cambia aquí para probar
H = 0.45       # altura/espesor (m)
E = 30e9       # módulo de Young (Pa)
rho = 2500.0
q_kNpm = 15.0  # carga distribuida (kN/m)
q = q_kNpm * 1e3  # N/m

# Discretización
N     = 80
Wdiv  = 8
Hdiv  = 4
SENSOR_X_FRACTION = 0.65

# Escalas visuales
DEF_SCALE      = 30.0
THICK_SCALE    = 1.0

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

        for i in range(N + 1):
            x = L * i / N
            ops.node(i + 1, x, 0.0)

        ops.fix(1,     1, 1, 0)
        ops.fix(N + 1, 1, 1, 0)

        A = B * H
        I = (B * (H**3)) / 12.0

        ops.geomTransf('Linear', 1)

        for i in range(1, N + 1):
            ops.element('elasticBeamColumn', i, i, i + 1, A, E, I, 1)

        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        p = -q * (L / N)
        for i in range(2, N):
            ops.load(i, 0.0, p, 0.0)

        ops.system('BandGeneral')
        ops.numberer('RCM')
        ops.constraints('Plain')
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        ops.analyze(1)

        w_def = [ops.nodeDisp(i + 1, 2) for i in range(N + 1)]
        return w_def

    except Exception:
        return None

def build_viz():
    w_def = deflection_opensees()
    if w_def is None:
        w_def = [deflection_analytic(L * i / N) for i in range(N + 1)]

    Ny = Wdiv + 1
    Nz = Hdiv + 1
    stride_x = Ny * Nz

    def vidx(i, j, k): return i * stride_x + j * Nz + k

    vertices, u_mag, indices = [], [], []

    halfB = B * 0.5
    halfH = (H * THICK_SCALE) * 0.5

    for i in range(N + 1):
        x = L * i / N
        w = w_def[i] * DEF_SCALE
        mm = abs(w_def[i]) * 1000.0
        for j in range(Ny):
            y = -halfB + (B * j / Wdiv)
            for k in range(Nz):
                h = -halfH + (2 * halfH * k / Hdiv)
                z = w + h
                vertices.extend([x, y, z])
                u_mag.append(mm)

    def add_quad(a, b, c, d): indices.extend([a, b, c, a, c, d])

    for i in range(N):
        for j in range(Wdiv):
            add_quad(vidx(i,j,Hdiv), vidx(i+1,j,Hdiv), vidx(i+1,j+1,Hdiv), vidx(i,j+1,Hdiv))
            add_quad(vidx(i,j,0), vidx(i+1,j,0), vidx(i+1,j+1,0), vidx(i,j+1,0))

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










