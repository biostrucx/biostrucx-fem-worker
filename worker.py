import os, time, math, json, traceback
from datetime import datetime, timezone
from pymongo import MongoClient, UpdateOne

# ====== Config desde ENV ======
MONGODB_URI   = os.getenv("MONGODB_URI")
MONGODB_DB    = os.getenv("MONGODB_DB", "biostrucx")
FEM_CLIENTS   = [c.strip() for c in os.getenv("FEM_CLIENTS", "jeimie").split(",") if c.strip()]
INTERVAL_SEC  = int(os.getenv("FEM_INTERVAL_SEC", "120"))  # cada 2 min por default

# Geometría/material por defecto (puedes parametrizar por ENV si quieres)
L = float(os.getenv("FEM_L", "25.0"))   # m
B = float(os.getenv("FEM_B", "1.0"))    # m
H = float(os.getenv("FEM_H", "1.0"))    # m
E = float(os.getenv("FEM_E", "30e9"))   # Pa
q = float(os.getenv("FEM_Q", "1000"))   # N/m (carga distribuida demo)

# Malla (resolución del ribbon)
N = int(os.getenv("FEM_N", "40"))       # nº segmentos a lo largo

# ====== Mongo ======
client = MongoClient(MONGODB_URI, serverSelectionTimeoutMS=15000)
db = client[MONGODB_DB]
col_result = db["simulation_result"]
col_ts     = db["simulation_ts"]

def ts_iso(dt=None):
    return (dt or datetime.now(timezone.utc)).isoformat()

def make_ribbon_vertices(L, B, H, u_func):
    """Crea un 'ribbon' (tira) para visualizar la viga y colorear por desplazamiento.
       Devuelve (vertices, indices, u_mag)."""
    verts = []
    indices = []
    u_mag = []
    halfW = 0.5 * B

    # malla 2 x (N+1) vértices (dos "bordes" del ribbon)
    for i in range(N + 1):
        x = L * i / N
        u = u_func(x)            # desplazamiento (m)
        y = 0.0
        z = 0.0
        # borde inferior/superior del ribbon
        verts.append([x, y, z])
        verts.append([x, y, z + halfW])
        u_mag.append(u)
        u_mag.append(u)

    # triángulos (dos por cada “cuadrado” a lo largo)
    for i in range(N):
        a = 2*i
        b = 2*i + 1
        c = 2*(i+1)
        d = 2*(i+1) + 1
        indices += [[a,b,c],[b,d,c]]

    return verts, indices, u_mag

def fem_beam_closed_form(L, E, B, H, q):
    """FEM 'mínimo viable': analítico de viga simplemente apoyada con carga distribuida."""
    I = (B * (H**3)) / 12.0
    def u(x):
        return (q * x * (L**3 - 2*L*(x**2) + x**3)) / (24.0 * E * I)
    u_mid = (5.0 * q * (L**4)) / (384.0 * E * I)  # deflexión máx (centro)
    return u, u_mid

def run_for_client(cid):
    # 1) Resuelve (analítico; si luego migras a OpenSeesPy, cambia aquí)
    u_func, u_mid = fem_beam_closed_form(L, E, B, H, q)

    # 2) Construye malla de visualización (ribbon) y colores por desplazamiento
    vertices, indices, u_mag = make_ribbon_vertices(L, B, H, u_func)

    # 3) Guarda "latest" para el 3D (simulation_result)
    now = datetime.now(timezone.utc)
    doc_result = {
        "clientid": cid,
        "ts": now,
        "status": "done",
        "params": {"desc":"viga simplemente apoyada", "L":L,"B":B,"H":H,"E":E,"q":q},
        "model": {"type":"beam"},
        "viz": {
            "vertices": vertices,   # [[x,y,z], ...]
            "indices":  indices,    # [[i,j,k], ...] triángulos
            "u_mag":    u_mag,      # color por desplazamiento
            "marker":   [L/2, 0, 0]  # posición del sensor (demo)
        }
    }
    col_result.update_one({"clientid": cid}, {"$set": doc_result}, upsert=True)

    # 4) Inserta punto en la serie (simulation_ts)
    col_ts.insert_one({
        "clientid": cid,
        "ts": now,
        "fem_mm": float(u_mid*1000.0)  # milímetros
    })

def main_loop():
    print("[worker] starting… clients=", FEM_CLIENTS, "every", INTERVAL_SEC, "sec")
    # conexión sanity-check
    client.admin.command("ping")
    while True:
        try:
            for cid in FEM_CLIENTS:
                run_for_client(cid)
        except Exception as e:
            traceback.print_exc()
        time.sleep(INTERVAL_SEC)

if __name__ == "__main__":
    main_loop()

