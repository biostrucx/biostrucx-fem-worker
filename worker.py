# --- construir una malla simple de viga (ribbon) y su deformada ---
import math
from datetime import datetime
from pymongo import MongoClient
import os

MONGODB_URI = os.environ["MONGODB_URI"]
MONGODB_DB  = os.environ.get("MONGODB_DB", "biostrucx")
CLIENTS     = os.environ.get("FEM_CLIENTS", "jeimie").split(",")

# parámetros básicos de demo (puedes ponerlos en env si quieres)
L   = 25.0   # largo [m]
W   = 1.0    # ancho visual del ribbon [m]
E   = 30e9   # [Pa]
I   = 0.25**4 / 12.0  # inercia aproximada (m^4) (25 cm ^4 /12, cambia a tu sección)
q   = 15000.0 # carga distribuida [N/m]
N   = 60     # divisiones a lo largo (malla Nx*1)

def beam_deflection_uniform(x, L, q, E, I):
    # viga simplemente apoyada con carga distribuida: w(x) = q x (L^3 - 2 L x^2 + x^3)/(24 E I)
    return q * x * (L**3 - 2*L*(x**2) + x**3) / (24.0 * E * I)

def make_beam_mesh(L=25.0, W=1.0, N=60):
    # malla rectangular de 2 filas (y = -W/2 y +W/2) y N segmentos a lo largo (x)
    xs = [L * i / N for i in range(N+1)]
    ys = [-W/2.0, +W/2.0]

    # vertices: (N+1)*2
    verts = []
    for x in xs:
        for y in ys:
            verts.extend([float(x), float(y), 0.0])  # z=0 plano

    # indices: 2 triángulos por cada quad a lo largo
    idx = []
    # índice (i fila y, j columna x):
    # v(i,j) = j*2 + i  (porque hay 2 filas en y)
    for j in range(N):
        v00 = j*2 + 0
        v01 = j*2 + 1
        v10 = (j+1)*2 + 0
        v11 = (j+1)*2 + 1
        # triángulos (v00, v10, v11) y (v00, v11, v01)
        idx.extend([v00, v10, v11,  v00, v11, v01])

    # u_mag por vértice (usa sólo x)
    u = []
    for j in range(N+1):
        x = xs[j]
        w = beam_deflection_uniform(x, L, q, E, I)  # [m]
        w_mm = 1000.0 * w                            # pasa a mm
        # dos filas en y: el mismo valor para ambos
        u.append(float(w_mm))  # para (j, fila 0)
        u.append(float(w_mm))  # para (j, fila 1)

    # validaciones útiles
    assert len(verts) % 3 == 0, "vertices debe ser múltiplo de 3"
    assert len(idx)   % 3 == 0, "indices debe ser múltiplo de 3"
    assert len(u) == len(verts)//3, "u_mag debe tener N_vertices"

    return verts, idx, u

def write_result_for(client_id: str):
    verts, idx, u = make_beam_mesh(L=L, W=W, N=N)

    doc = {
        "clientid": client_id,
        "status":   "done",
        "ts":       datetime.utcnow(),
        "model":    {"type": "beam_ribbon", "note": "demo uniform load"},
        "params":   {"L": L, "W": W, "E": E, "I": I, "q": q, "N": N},
        "viz": {
            "vertices": verts,        # lista plana [x,y,z,...]
            "indices":  idx,          # lista plana [i0,i1,i2,...]
            "u_mag":    u,            # mm para color
            "marker":   [L*0.25, 0.0, 0.0]  # donde está el sensor (demo)
        }
    }

    cli = MongoClient(MONGODB_URI)
    db  = cli[MONGODB_DB]
    db.simulation_result.insert_one(doc)

    # también escribe un punto a la serie temporal (para el gráfico 1)
    db.simulation_ts.insert_one({
        "clientid": client_id,
        "ts": datetime.utcnow(),
        "fem_mm": float(max(u))  # por ejemplo el máximo en el paso
    })

    print(f"[worker] wrote viz: V={len(verts)//3} verts, T={len(idx)//3} tris, max={max(u):.3f} mm")

# --- bucle simple (o deja tu cron actual) ---
if __name__ == "__main__":
    for cid in [c.strip() for c in CLIENTS if c.strip()]:
        write_result_for(cid)

