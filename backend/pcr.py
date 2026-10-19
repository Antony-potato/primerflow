"""PCR virtual: pares contra secuencias, ambas hebras, con desajustes (docs/MOTOR_CALCULO.md §9).

Una sola regla cubre pares normales, invertidos, de un solo primer y cruzados del multiplex:
hay producto si el primer x se pega en s (sitio directo, en i) y el primer y se pega como rc(y)
más adelante (termina en j), con j - i <= tamano_max."""
from collections import Counter

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

import motor


def _arr(s):
    return np.frombuffer(s.encode(), dtype=np.uint8)


def _sitios(s, p, max_mm, estricto, extremo3_al_final):
    k = len(p)
    if len(s) < k:
        return np.empty(0, int), np.empty(0, int)
    dif = sliding_window_view(s, k) != p
    mm = dif.sum(1)
    ok = mm <= max_mm
    if estricto:  # 0 desajustes en los 3 nt del extremo 3′
        ok &= ~(dif[:, -3:] if extremo3_al_final else dif[:, :3]).any(1)
    idx = np.flatnonzero(ok)
    return idx, mm[idx]


def productos(primers, seq, max_mm=2, tamano_max=3000, estricto=False):
    """primers: [(nombre, seq)]. Devuelve [(x, y, inicio, tamaño, desajustes_x, desajustes_y)]."""
    s = _arr(seq.upper().replace("-", ""))
    directos = {n: _sitios(s, _arr(p), max_mm, estricto, True) for n, p in primers}
    inversos = {n: _sitios(s, _arr(motor.rc(p)), max_mm, estricto, False) for n, p in primers}
    largo = dict((n, len(p)) for n, p in primers)
    out = []
    for nx, (ix, mx) in directos.items():
        for ny, (iy, my) in inversos.items():
            fin = iy + largo[ny]
            for i, a in zip(ix, mx):
                ok = (iy >= i) & (fin - i <= tamano_max)
                out += [(nx, ny, int(i), int(j - i), int(a), int(b)) for j, b in zip(fin[ok], my[ok])]
    return out


def pcr_virtual(pares, secuencias, max_mm=2, tamano_max=3000, estricto=False, cruzados=True):
    resultados, resumen = [], []
    for par in pares:
        nF, nR = par["nombre"] + "_cF", par["nombre"] + "_cR"
        primers = [(nF, motor.limpiar(par["F"])), (nR, motor.limpiar(par["R"]))]
        propias = [0, 0]
        otras, tamanos = 0, Counter()
        for sq in secuencias:
            prods = productos(primers, sq["seq"], max_mm, tamano_max, estricto)
            mejor = min(prods, key=lambda p: (p[4] + p[5], p[0] != nF)) if prods else None
            mismo = sq["especie"] == par["especie"]
            if mismo:
                propias[1] += 1
                propias[0] += bool(prods)
            else:
                otras += bool(prods)
            if mejor:
                tamanos[mejor[3]] += 1
            resultados.append({
                "par": par["nombre"], "secuencia": sq["id"], "especie": sq["especie"], "esperado": mismo,
                "producto": bool(prods), "tamano": mejor[3] if mejor else None,
                "desajustes": f"{mejor[4]}/{mejor[5]}" if mejor else None,
                "combinacion": f"{mejor[0]} + {mejor[1]}" if mejor else None, "otros_productos": max(0, len(prods) - 1),
            })
        resumen.append({"par": par["nombre"], "especie": par["especie"], "tubo": par.get("tubo"),
                        "propia": f"{propias[0]}/{propias[1]}", "otras": otras,
                        "probadas_otras": len(secuencias) - propias[1], "tamanos": sorted(tamanos)})

    cruces = []
    if cruzados:  # F de un par con R (o F) de otro, dentro del mismo tubo
        tubos = {}
        for p in pares:
            tubos.setdefault(p.get("tubo") or "", []).append(p)
        for tubo, ps in tubos.items():
            dueno = {}
            pool = []
            for p in ps:
                for lado in ("F", "R"):
                    nombre = f"{p['nombre']}_c{lado}"
                    dueno[nombre] = p["nombre"]
                    pool.append((nombre, motor.limpiar(p[lado])))
            for sq in secuencias:
                for x, y, i, t, a, b in productos(pool, sq["seq"], max_mm, tamano_max, estricto):
                    if dueno[x] != dueno[y]:
                        cruces.append({"tubo": tubo, "secuencia": sq["id"], "especie": sq["especie"],
                                       "primers": f"{x} + {y}", "tamano": t, "desajustes": f"{a}/{b}"})
    return {"resultados": resultados, "resumen": resumen, "cruzados": cruces}
