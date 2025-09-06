#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BioStrucX – FEM Worker (demo beam)
- Calcula una viga simplemente apoyada con carga distribuida (analítico).
- Guarda:
    * biostrucx.simulation_result: malla + colores (viz) y metadatos
    * biostrucx.simulation_ts: serie temporal fem_mm

Variables de entorno:
    MONGODB_URI           (obligatoria)
    MONGODB_DB            (por defecto: biostrucx)
    FEM_CLIENTS           (coma-sep, ej: "jeimie,client2"; por defecto: "jeimie")
    FEM_INTERVAL_SEC      (por defecto: 120)

Parámetros (opcionales) del modelo:
    FEM_L, FEM_B, FEM_H, FEM_E, FEM_Q, FEM_N, FEM_W
"""

import os
import time
import math
from datetime import datetime, timezone

from pymongo import MongoClient, ASCENDING, DESCENDING
from pymongo.errors import OperationFailure, PyMongoError

# ==========================
# Config de entorno
# ==========================
MONGO_URI = os.getenv("MONGODB_URI") or os.getenv("MONGODB_URI".lower()) or os.getenv("mongodb_uri")
MONGO_DB = os.getenv("MONGODB_DB", "biostrucx")

CLIENTS = [c.strip() for c in (os.getenv("FEM_CLIENTS") or "jeimie").split(",") if c.strip()]
INTERVAL_SEC = int(os.getenv("FEM_INTERVAL_SEC", "120"))

# Parámetros de la viga (unidades SI)
L = float(os.getenv("FEM_L", "25.0"))           # longitud (m)
B = float(os.getenv("FEM_B", "0.25"))           # base sección (m)
H = float(os.getenv("FEM_H", "0.25"))           # altura sección (m)
E = float(os.getenv("FEM_E", "2.1e11"))         # módulo E (Pa)
Q = float(os.getenv("FEM_Q", "15000"))          # carga distribuida (N/m), p.ej. 15 kN/m
N = int(os.getenv("FEM_N", "80"))               # segmentos de malla (>= 2)
W = float(os.getenv("FEM_W", "0.40"))           # ancho visual del “ribbon” (m) para render

# ==========================
# Conexión Mongo + índices
# ==========================
def get_db():
    if not MONGO_URI:
        raise RuntimeError("Falta MONGODB_URI")

    cli = MongoClient(MONGO_URI)
    db = cli[MONGO_DB]

    # Crea el MISMO nombre de índice en ambas colecciones.
    # Si ya existe con otro nombre, capturamos code=85 (IndexOptionsConflict) y seguimos.
    try:
        db.simulation_result.create_index(
            [("clientid", ASCENDING), ("ts", DESCENDING)],
            name="idx_clientid_ts",
            background=True,
        )
    except OperationFailure as e:
        if getattr(e, "code", None) != 85:
            raise

    try:
        db.simulation_ts.create_index(
            [("clientid", ASCENDING), ("ts", DESCENDING)],
            name="idx_clientid_ts",
            background=True,
        )
    except OperationFailure as e:
        if getattr(e, "code", None) != 85:
            raise

    return db


# ==========================
# Modelo demo de viga (analítico)
# ==========================
def beam_deflection_udl(x, L, E, I, q):
    """
    Deflexión w(x) de viga simplemente apoyada con carga distribuida uniforme q (N/m).
    E en Pa, I en m^4, resultado en metros.
    Fórmula clásica: w(x) = q x (L^3 - 2 L x^2 + x^3) / (24 E I)
    """
    return (q * x * (L**3 - 2.0 * L * (x**2) + x**3)) / (24.0 * E * I)


def build_beam_ribbon(L, W, N, E, B, H, q):
    """
    Construye una “malla cinta” (ribbon) a lo largo de X, ancho W (en Y).
    - vertices: [x,y,z, x,y,z, ...]  (len % 3 == 0)
    - indices: [i0,i1,i2, ...]       (triángulos; len % 3 == 0)
    - u_mag:   [mm por vértice]      (len == num_vertices/3)
    - marker:  [xm, ym, zm]          (punto de interés; aquí, midspan)
    La coordenada z se deforma con la deflexión analítica.
    """
    N = max(2, int(N))
    I = (B * (H**3)) / 12.0  # inercia rectangular

    xs = [L * i / N for i in range(N + 1)]
    half_w = W / 2.0

    vertices = []
    u_mag = []

    # 2 vértices por "sección": y = -half_w y y = +half_w
    for x in xs:
        w = beam_deflection_udl(x, L, E, I, q)  # metros (positivo hacia abajo)
        z = -w                                  # sign: hacia abajo negativo (opcional)
        mm = w * 1000.0

        # lado -W/2
        vertices.extend([x, -half_w, z])
        u_mag.append(float(mm))

        # lado +W/2
        vertices.extend([x, +half_w, z])
        u_mag.append(float(mm))

    # Triángulos entre dos franjas consecutivas
    # Para cada i: (i*2, i*2+1, i*2+2) y (i*2+1, i*2+3, i*2+2)
    indices = []
    for i in range(N):
        a = i * 2
        b = a + 1
        c = a + 2
        d = a + 3
        indices.extend([a, b, c])
        indices.extend([b, d, c])

    # Punto de interés (midspan)
    xmid = L / 2.0
    wmid = beam_deflection_udl(xmid, L, E, I, q)
    marker = [xmid, 0.0, -wmid]

    # Validaciones mínimas (evitan subir algo “roto”)
    assert len(vertices) % 3 == 0, "vertices debe ser múltiplo de 3"
    assert len(indices) % 3 == 0, "indices debe ser múltiplo de 3"
    assert len(u_mag) == len(vertices) // 3, "u_mag debe tener N_vértices"

    # Desplazamiento (mm) en el medio de la luz, útil para la serie temporal
    mid_mm = float(wmid * 1000.0)
    return vertices, indices, u_mag, marker, mid_mm


# ==========================
# Ejecución por cliente
# ==========================
def run_for_client(db, clientid: str):
    """
    Calcula la viga demo, guarda en:
      - simulation_result (documento con viz)
      - simulation_ts (punto fem_mm)
    """
    ts_now = datetime.now(timezone.utc)

    try:
        vertices, indices, u_mag, marker, mid_mm = build_beam_ribbon(
            L=L, W=W, N=N, E=E, B=B, H=H, q=Q
        )

        doc_result = {
            "clientid": clientid,
            "model": {"type": "beam_demo"},
            "params": {
                "L": L, "B": B, "H": H, "E": E, "q": Q, "N": N, "W": W
            },
            "status": "done",
            "ts": ts_now,
            "viz": {
                "vertices": [float(v) for v in vertices],
                "indices": [int(i) for i in indices],
                "u_mag":   [float(v) for v in u_mag],
                "marker":  [float(x) for x in marker],
            },
        }

        # Insert resultado + timeseries
        db.simulation_result.insert_one(doc_result)
        db.simulation_ts.insert_one({
            "ts": ts_now,
            "clientid": clientid,
            "fem_mm": float(mid_mm),
        })

        print(f"[worker] {clientid} ok | fem_mm={mid_mm:.3f} mm")

    except AssertionError as e:
        print(f"[worker] {clientid} viz assertion error: {e}")
        db.simulation_result.insert_one({
            "clientid": clientid,
            "status": "error",
            "error": f"viz_assert: {str(e)}",
            "ts": ts_now,
        })
    except PyMongoError as e:
        print(f"[worker] {clientid} mongo error: {e}")
    except Exception as e:
        print(f"[worker] {clientid} general error: {e}")
        db.simulation_result.insert_one({
            "clientid": clientid,
            "status": "error",
            "error": f"runtime: {str(e)}",
            "ts": ts_now,
        })


# ==========================
# Main loop
# ==========================
def main():
    print(f"[worker] starting.. clients={CLIENTS} every {INTERVAL_SEC} sec")
    db = get_db()

    while True:
        start = time.time()
        for cid in CLIENTS:
            run_for_client(db, cid)

        # Sleep hasta completar el intervalo, respetando tiempo de ejecución
        spent = time.time() - start
        to_sleep = max(1.0, INTERVAL_SEC - spent)
        time.sleep(to_sleep)


if __name__ == "__main__":
    main()
