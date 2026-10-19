"""Alineamiento, conservación (estilo MEGA / MultAlin), diseño con primer3 y ranking.
Coordenadas internas: base 0, intervalos [inicio, fin)."""
import io
import re
from collections import Counter

import primer3
from Bio import AlignIO, SeqIO

import motor
import pcr

LETRAS = set("ACGTNRYKMBVDHSW-")


def leer(texto, formato="fasta"):
    """FASTA, o el .mul (FASTA alineado) / .msf de MultAlin."""
    if formato == "msf":
        return [{"id": r.id, "seq": str(r.seq)} for r in AlignIO.read(io.StringIO(texto), "msf")]
    return [{"id": r.id, "seq": str(r.seq)} for r in SeqIO.parse(io.StringIO(texto), "fasta")]


def _msa(seqs):
    try:
        from pyfamsa import Aligner, Sequence
    except ImportError:
        raise ValueError("Las secuencias miden distinto y falta pyfamsa (pip install pyfamsa). "
                         "También puedes subir el alineamiento de MultAlin (.mul o .msf).")
    msa = Aligner(guide_tree="upgma").align([Sequence(str(i).encode(), s.replace("-", "").encode())
                                             for i, s in enumerate(seqs)])
    por_id = {int(s.id.decode()): s.sequence.decode() for s in msa}
    return [por_id[i] for i in range(len(seqs))]


def alinear(secuencias, ya_alineado=False, alto=1.0, bajo=0.5):
    ids = [s["id"] for s in secuencias]
    seqs = [re.sub(r"\s", "", s["seq"]).upper().replace(".", "-").replace("~", "-") for s in secuencias]
    if len(seqs) < 2:
        raise ValueError("Se necesitan al menos 2 secuencias")
    for i, s in zip(ids, seqs):
        if not s or not set(s) <= LETRAS:
            raise ValueError(f"{i}: caracteres no válidos")
    if len(set(map(len, seqs))) > 1:
        if ya_alineado:
            raise ValueError("El alineamiento importado tiene filas de distinto largo")
        filas = _msa(seqs)
    else:
        filas = seqs  # ponytail: mismo largo = ya alineadas (válido sin indels; si no, subir el de MultAlin)

    n = len(filas)
    consenso, multalin, cons_col, variables = [], [], [], []
    C = V = Pi = hueco = 0
    for col in range(len(filas[0])):
        letras = [f[col] for f in filas]
        cnt = Counter(letras)
        base, k = max(((b, c) for b, c in cnt.items() if b != "-"), key=lambda x: x[1], default=("-", 0))
        if k < cnt.get("-", 0):
            base = "-"
        f = k / n if base != "-" else 0.0
        consenso.append(base)
        cons_col.append(round(f, 3))
        multalin.append(base if f >= alto else base.lower() if f >= bajo else ".")
        if "-" in cnt:
            hueco += 1
        elif len(cnt) == 1:
            C += 1
        else:
            V += 1
            Pi += sum(c >= 2 for c in cnt.values()) >= 2
        if len(cnt) > 1:
            variables.append({"columna": col, "mayoritaria": base, "composicion": dict(cnt.most_common()),
                              "cepas": {ids[r]: x for r, x in enumerate(letras) if x != base}})

    mapa = [c for c, b in enumerate(consenso) if b != "-"]  # base de la plantilla → columna
    conservacion = [cons_col[c] for c in mapa]
    bloques, ini = [], None
    for i, f in enumerate(conservacion + [0]):
        if f == 1 and ini is None:
            ini = i
        elif f < 1 and ini is not None:
            if i - ini >= 18:
                bloques.append({"inicio": ini, "fin": i, "largo": i - ini})
            ini = None
    bloques.sort(key=lambda b: -b["largo"])
    return {
        "filas": [{"id": i, "seq": f} for i, f in zip(ids, filas)],
        "consenso": "".join(consenso), "consenso_multalin": "".join(multalin), "conservacion_columna": cons_col,
        "plantilla": "".join(consenso[c] for c in mapa), "conservacion": conservacion, "mapa": mapa,
        "variables": variables, "bloques": bloques,
        "resumen": {"secuencias": n, "total": len(consenso), "con_hueco": hueco, "analizables": len(consenso) - hueco,
                    "C": C, "V": V, "Pi": Pi, "S": V - Pi},
    }


def _restar(rangos, excluir):
    out = []
    for a, b in rangos:
        trozos = [(a, b)]
        for c, d in excluir:
            nuevos = []
            for x, y in trozos:
                if d < x or c > y:
                    nuevos.append((x, y))
                    continue
                if x < c:
                    nuevos.append((x, c - 1))
                if d < y:
                    nuevos.append((d + 1, y))
            trozos = nuevos
        out += [[x, y] for x, y in trozos if y - x >= 10]
    return out


def _cruzados(F, R, oligos_tubo, plantillas):
    """Productos del candidato con primers de otros pares del tubo (PCR virtual sobre cada plantilla del tubo)."""
    pool = [("F", F), ("R", R)] + [(o["nombre"], motor.limpiar(o["seq"])) for o in oligos_tubo]
    propios = {"F", "R"}
    return [f"{x} + {y} en {nombre} ({t} pb)" for nombre, seq in plantillas
            for x, y, _, t, _, _ in pcr.productos(pool, seq) if (x in propios) != (y in propios)]


def evaluar_par(prefijo, f_seq, r_seq, f_sitio, r_sitio, plantilla, conservacion, cond, oligos_tubo=(), plantillas_tubo=()):
    F = motor.informe(prefijo + "_cF", f_seq, cond)
    R = motor.informe(prefijo + "_cR", r_seq, cond)
    F.update(inicio=f_sitio[0], fin=f_sitio[1], hebra=1)
    R.update(inicio=r_sitio[0], fin=r_sitio[1], hebra=-1)
    pcr = cond["pcr"]
    hetero = motor.dimero(F["seq"], R["seq"], pcr)

    var_f = [p for p in range(*f_sitio) if conservacion[p] < 1]
    var_r = [p for p in range(*r_sitio) if conservacion[p] < 1]
    en_3p = [p for p in var_f if p >= f_sitio[1] - 5] + [p for p in var_r if p < r_sitio[0] + 5]

    tubo = [(o["nombre"], motor.dimero(x, motor.limpiar(o["seq"]), cond.get("multiplex", pcr))["dg"])
            for o in oligos_tubo for x in (F["seq"], R["seq"])]
    peor = min(tubo, key=lambda t: t[1], default=(None, None))
    cruces = _cruzados(F["seq"], R["seq"], oligos_tubo, [(prefijo, plantilla)] + list(plantillas_tubo)) if oligos_tubo else []

    par = {
        "cruzados": len(cruces), "cruzados_detalle": cruces[:10],
        "F": F, "R": R, "dtm": round(abs(F["tm_pcr"] - R["tm_pcr"]), 1),
        "hetero_dg": hetero["dg"], "anclaje_3p": motor.anclaje_3p(F["seq"], R["seq"], pcr),
        "tubo_dg": peor[1], "tubo_con": peor[0], "variables": len(var_f) + len(var_r), "variables_3p": len(en_3p),
        "estructuras": {"heterodimero": hetero},
        "producto": {"inicio": f_sitio[0], "fin": r_sitio[1], "tamano": r_sitio[1] - f_sitio[0],
                     "gc_pct": round(motor.gc_pct(plantilla[f_sitio[0]:r_sitio[1]]), 1),
                     "seq": plantilla[f_sitio[0]:r_sitio[1]]},
    }
    par["niveles"], par["banderas"] = motor.evaluar(par, motor.CRITERIOS_PAR)
    if en_3p:
        par["niveles"]["variables"] = "rojo"
        par["banderas"] = [b for b in par["banderas"] if b["campo"] != "variables"] + [
            {"regla": "Posición variable en los 5 nt 3′", "campo": "variables", "valor": len(en_3p), "nivel": "rojo"}]
    todas = [("F · " + b["regla"], b) for b in F["banderas"]] + [("R · " + b["regla"], b) for b in R["banderas"]] + \
            [(b["regla"], b) for b in par["banderas"]]
    par["rojos"] = [f"{k}: {b['valor']}" for k, b in todas if b["nivel"] == "rojo"]
    amarillos = sum(b["nivel"] == "amarillo" for _, b in todas)

    # ponytail: pesos heurísticos (LOGICA_PRODUCTO.md §5); pasar al perfil cuando haya que ajustarlos
    d = [
        ("Tm F lejos de 60.5 °C", 3 * abs(F["tm_pcr"] - 60.5)),
        ("Tm R lejos de 60.5 °C", 3 * abs(R["tm_pcr"] - 60.5)),
        ("ΔTm F–R", 5 * par["dtm"]),
        ("%GC lejos de 50", 0.3 * (abs(F["gc_pct"] - 50) + abs(R["gc_pct"] - 50))),
        ("Heterodímero F×R", 2 * max(0, -par["hetero_dg"] - 3)),
        ("Autodímeros", max(0, -F["autodimero_dg"] - 3) + max(0, -R["autodimero_dg"] - 3)),
        ("Horquillas", max(0, -F["horquilla_dg"] - 1) + max(0, -R["horquilla_dg"] - 1)),
        ("Dímeros con el tubo", 2 * max(0, -(peor[1] or 0) - 4)),
        ("Avisos en amarillo", 8 * amarillos),
    ]
    par["desglose"] = [{"regla": k, "puntos": -round(p, 1)} for k, p in d if p >= 0.05]
    par["puntaje"] = max(0, round(100 - sum(p for _, p in d)))
    return par


def disenar(plantilla, conservacion, prefijo, cond, region=None, tamanos=((150, 600),), excluir=(),
            oligos_tubo=(), n=20, min_conservacion=100, plantillas_tubo=()):
    if len(plantilla) != len(conservacion):
        raise ValueError("plantilla y conservación deben medir lo mismo")
    rangos = _restar(tamanos, excluir)
    if not rangos:
        raise ValueError("No queda ningún rango de tamaño libre en el tubo")
    pcr = cond["pcr"]
    seq_args = {"SEQUENCE_ID": prefijo, "SEQUENCE_TEMPLATE": plantilla,
                "SEQUENCE_QUALITY": [round(100 * c) for c in conservacion]}
    if region:
        seq_args["SEQUENCE_INCLUDED_REGION"] = [region[0], region[1] - region[0]]
    r = primer3.design_primers(seq_args, {
        "PRIMER_NUM_RETURN": 50, "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 24,
        "PRIMER_OPT_TM": 60.5, "PRIMER_MIN_TM": 57.0, "PRIMER_MAX_TM": 64.0,
        "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
        "PRIMER_GC_CLAMP": 0, "PRIMER_MAX_END_GC": 3, "PRIMER_MAX_POLY_X": 4,
        "PRIMER_PAIR_MAX_DIFF_TM": 2.0, "PRIMER_MAX_END_STABILITY": 9.0,
        "PRIMER_SALT_MONOVALENT": pcr["na_mM"], "PRIMER_SALT_DIVALENT": pcr["mg_mM"],
        "PRIMER_DNTP_CONC": pcr["dntp_mM"], "PRIMER_DNA_CONC": pcr["oligo_uM"] * 1000,
        "PRIMER_TM_FORMULA": 1, "PRIMER_SALT_CORRECTIONS": 2, "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_MIN_QUALITY": min_conservacion, "PRIMER_MIN_END_QUALITY": 100,
        "PRIMER_QUALITY_RANGE_MIN": 0, "PRIMER_QUALITY_RANGE_MAX": 100,
        "PRIMER_PRODUCT_SIZE_RANGE": rangos,
    })
    candidatos, descartados = [], []
    for i in range(r.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        fp, fl = r[f"PRIMER_LEFT_{i}"]
        rp, rl = r[f"PRIMER_RIGHT_{i}"]  # primer3 da la posición del extremo 5′ del reverse (la más a la derecha)
        par = evaluar_par(prefijo, r[f"PRIMER_LEFT_{i}_SEQUENCE"], r[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                          (fp, fp + fl), (rp - rl + 1, rp + 1), plantilla, conservacion, cond, oligos_tubo, plantillas_tubo)
        par["id"] = f"{prefijo}-{i + 1}"
        par["penalizacion_primer3"] = round(r[f"PRIMER_PAIR_{i}_PENALTY"], 2)
        (descartados if par["rojos"] else candidatos).append(par)
    candidatos.sort(key=lambda p: (-p["puntaje"], p["penalizacion_primer3"]))
    explicacion = {k: r[k] for k in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN") if k in r}
    return {"candidatos": candidatos[:n], "descartados": descartados[:20], "rangos": rangos, "explicacion": explicacion}
