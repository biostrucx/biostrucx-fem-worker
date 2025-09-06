import os
import time
import math
from datetime import datetime, timezone
from typing import Callable, List, Tuple

from pymongo import MongoClient, ASCENDING, DESCENDING
from pymongo.errors import PyMongoError

# --- Utilidades de entorno ----------------------------------------------------

def getenv_str(key: str, default: str) -> str:
    v = os.getenv(key)
    return v if (v is not None and v != "") else default

def getenv_float(key: str, default: float) -> float:
    v = os.getenv(key)
    try:
        return float(v) if v not in (None, "") else default
    except Exception:
        return default

def getenv_int(key: str, default: int) -> int:
    v = os.getenv(key)
    try:
        return int(v) if v not in (None, "") else default
    except Exception:
        return default

# --- Parámetros por defecto (puedes sobreescribir con env) --------------------

FEM_L      = getenv_float("FEM_L",      25.0)        # Longitud (m)
FEM_B      = getenv_float("FEM_B",       0.25)        # base sección (m)
FEM_H      = getenv_float("FEM_H",       0.25)        # altura sección (m)
FEM_E      = getenv_float("FEM_E",  200e9)            # Módulo Young (Pa)
FEM_Q      = getenv_float("FEM_Q",   15e3)            # Carga distribuida (N/m)
FEM_N      = getenv_int  ("FEM_N",        80)         # Segmentos a lo largo
FEM_HALFW  = getenv_float("FEM_HALFW",   0.4)         # semiancho del ribbon (visual)
INTERVAL_S = getenv_int  ("FEM_INTERVAL_SEC", 120)

MONGO_URI  = getenv_str("MONGODB_URI", "")
MONGO_DB   = getenv_str("MONGODB_DB",  "biostrucx")
CLIENTS    = [c.strip() for c in getenv_str("FEM_CLIENTS", "jeimie").split(",") if c.strip()]

# --- Mongo --------------------------------------------------------------------

def get_db():
    if not MONGO_URI:
        raise RuntimeError("Falta MONGODB_URI")
    cli = MongoClient(MONGO_URI)
    db = cli[MONGO_DB]
    # Índices útiles
    db.simulation_result.create_index([("clientid", ASCENDING), ("ts", DESCENDING)])
    db.simulation_ts.create_index([("clientid", ASCENDING), ("ts", DESCENDING)])
    return db

# --- Analítica Euler-Bernoulli (viga simplemente apoyada, carga distribuida) --

def euler_bernoulli_udl(L: float, E: float, I: float, q: float) -> Tuple[Callable[[float], float], float]:
    """
    Devuelve:
      w(x): deflexión en metros para x∈[0,L]
      w_max: deflexión máxima en el centro (m)
    Fórmula clásica UDL en viga simplemente apoyada.
    """
    def w(x: float) -> float:
        # w(x) = (q x (L^3 - 2 L x^2 + x^3)) / (24 E I)
        return (q * x * (L**3 - 2*L*(x**2) + x**3)) / (24.0 * E * I)

    w_max = 5.0 * q * (L**4) / (384.0 * E * I)  # en el centro L/2
    return w, w_max

# --- Generación de malla tipo "ribbon" ----------------------------------------

def make_beam_ribbon(
    L: float,
    N: int,
    halfW: float,
    u_func: Callable[[float], float]
) -> Tuple[List[float], List[int], List[float]]:
    """
    Genera una “cinta” (2 vértices por estación a lo largo del eje X).
    Devuelve:
      vertices: [x,y,z, x,y,z, ...]  (longitud múltiplo de 3)
      indices:  [i0,i1,i2, i0,i1,i2, ...] (múltiplo de 3)
      u_mag:    [u0,u1, ...] por vértice
    """
    xs = [ (i * L) / N for i in range(N + 1) ]  # N+1 estaciones

    vertices: List[float] = []
    u_mag:   List[float] = []

    for x in xs:
        # Ribbon “plano”; sólo coloreamos por u_mag.
        # (Si quisieras deformar, podrías poner z = k * w(x))
        yL, yR, z = -halfW, +halfW, 0.0
        w = u_func(x)                       # metros
        u_mm = abs(w) * 1000.0              # para color, en mm

        # 2 vértices por estación
        vertices.extend([x, yL, z])
        u_mag.append(u_mm)

        vertices.extend([x, yR, z])
        u_mag.append(u_mm)

    # Triangulación por quad: 2 triángulos por segmento
    indices: List[int] = []
    for k in range(N):
        a = 2*k
        b = 2*k + 1
        c = 2*(k+1)
        d = 2*(k+1) + 1
        # tri 1
        indices.extend([a, b, c])
        # tri 2
        indices.extend([c, b, d])

    # Validaciones duras
    assert len(vertices) % 3 == 0, f"vertices len {len(vertices)} no múltiplo de 3"
    assert len(indices)  % 3 == 0, f"indices len {len(indices)} no múltiplo de 3"
    Nv = len(vertices) // 3
    assert len(u_mag) == Nv, f"u_mag {len(u_mag)} != Nv {Nv}"
    # Índices en rango
    if indices:
        mi, ma = min(indices), max(indices)
        assert 0 <= mi and ma < Nv, f"indices fuera de rango [0,{Nv-1}]"

    return vertices, indices, u_mag

# --- Worker principal ---------------------------------------------------------

def run_once_for_client(db, client_id: str):
    # Propiedades geométricas
    L = FEM_L
    B = FEM_B
    H = FEM_H
    E = FEM_E
    q = FEM_Q
    N = FEM_N
    halfW = FEM_HALFW

    I = (B * (H**3)) / 12.0   # inercia rectangular

    # Deflexión analítica
    w_func, w_max = euler_bernoulli_udl(L, E, I, q)

    # Malla “ribbon”
    vertices, indices, u_mag = make_beam_ribbon(L, N, halfW, w_func)

    Nv = len(vertices)//3
    Ntri = len(indices)//3
    print(f"[worker] client={client_id} Nv={Nv} ({len(vertices)} floats), Ntri={Ntri}, Nu={len(u_mag)}")

    # Marker (p.ej. centro de la viga)
    marker = [L/2.0, 0.0, 0.0]

    # Documento principal (para el visor)
    doc = {
        "clientid": client_id,
        "status": "done",
        "ts": datetime.now(timezone.utc),
        "model": {
            "type": "beam_demo",
            "note": "Euler-Bernoulli simply supported, UDL"
        },
        "params": {
            "L": L, "B": B, "H": H, "E": E, "q": q, "N": N
        },
        "viz": {
            "vertices": vertices,
            "indices": indices,
            "u_mag": u_mag,
            "marker": marker
        }
    }

    # Escribe resultado “snapshot”
    db.simulation_result.insert_one(doc)

    # Serie de tiempo (deflexión en el centro, en mm)
    fem_mm = abs(w_max) * 1000.0
    db.simulation_ts.insert_one({
        "clientid": client_id,
        "ts": datetime.now(timezone.utc),
        "fem_mm": float(fem_mm)
    })

def main():
    print("[worker] starting… clients:", CLIENTS, f"every {INTERVAL_S} sec")
    db = get_db()

    while True:
        for cid in CLIENTS:
            try:
                run_once_for_client(db, cid)
            except AssertionError as ae:
                # No subas docs rotos: loggea y sigue
                print(f"[worker][{cid}] ASSERT FAIL:", str(ae))
            except PyMongoError as me:
                print(f"[worker][{cid}] Mongo error:", str(me))
            except Exception as ex:
                print(f"[worker][{cid}] Error:", str(ex))
        try:
            time.sleep(INTERVAL_S)
        except KeyboardInterrupt:
            break

if __name__ == "__main__":
    main()

