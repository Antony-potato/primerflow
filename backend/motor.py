"""Motor de cálculo estilo IDT OligoAnalyzer. Fórmulas y validación: docs/MOTOR_CALCULO.md."""
import math
import re

import primer3
from Bio.SeqUtils import MeltingTemp as mt

VERSION = "0.1.0"

PRESETS = {
    "fabrica": {"oligo_uM": 0.25, "na_mM": 50, "mg_mM": 0, "dntp_mM": 0},
    "pcr": {"oligo_uM": 0.4, "na_mM": 50, "mg_mM": 2, "dntp_mM": 0.2},
    "multiplex": {"oligo_uM": 0.2, "na_mM": 50, "mg_mM": 2.5, "dntp_mM": 0.2},
}

MW = dict(A=313.209, C=289.184, G=329.208, T=304.196)  # IDT, residuo anhidro
EPS_NN = dict(AA=27400, AC=21200, AG=25000, AT=22800, CA=21200, CC=14600, CG=18000, CT=15200,
              GA=25200, GC=17600, GG=21600, GT=20000, TA=23400, TC=16200, TG=19000, TT=16800)
EPS_1 = dict(A=15400, C=7400, G=11500, T=8700)
COMP = str.maketrans("ACGTNRYKMBVDHSW-", "TGCANYRMKVBHDSW-")
INF = math.inf

# Perfil "Tesis" (docs/LOGICA_PRODUCTO.md §4): (etiqueta, campo, verde [min, max], amarillo [min, max]).
# Fuera del amarillo = rojo = rechazo.
CRITERIOS_OLIGO = [
    ("Largo (nt)", "largo", (18, 24), (16, 28)),
    ("%GC", "gc_pct", (40, 60), (35, 65)),
    ("Tm fábrica (°C)", "tm_fabrica", (50, 55), (48, 57)),
    ("Tm PCR (°C)", "tm_pcr", (58, 63), (57, 64)),
    ("Pinza GC (G/C en los 5 nt 3′)", "pinza_gc", (1, 2), (1, 3)),
    ("Repetición más larga (nt)", "poly_x", (0, 3), (0, 4)),
    ("Horquilla ΔG 25 °C", "horquilla_dg", (-3, INF), (-6, INF)),
    ("Horquilla Tm (°C)", "horquilla_tm", (-INF, 47), (-INF, 57)),
    ("Autodímero ΔG PCR", "autodimero_dg", (-5, INF), (-6, INF)),
    ("Estabilidad 3′ (ΔG 5 nt)", "estabilidad_3p", (-9, INF), (-9, INF)),
]
CRITERIOS_PAR = [
    ("ΔTm F–R (°C)", "dtm", (0, 1), (0, 2)),
    ("Heterodímero F×R ΔG", "hetero_dg", (-6, INF), (-9, INF)),
    ("Anclaje 3′ F×R", "anclaje_3p", (-5, INF), (-5, INF)),
    ("Peor dímero con el tubo", "tubo_dg", (-6, INF), (-9, INF)),
    ("Productos cruzados en el tubo", "cruzados", (0, 0), (0, 0)),
    ("Posiciones variables", "variables", (0, 0), (0, INF)),  # rojo si hay alguna en los 5 nt 3′ (diseno.py)
]


def rc(s):
    return s.translate(COMP)[::-1]


def limpiar(seq):
    s = re.sub(r"[\s\d]", "", seq).upper()
    if not s or not set(s) <= set("ACGT"):
        raise ValueError(f"Secuencia inválida (sólo A, C, G, T): {seq[:40]}")
    return s


def tm(s, c):
    # Ct/4 + Owczarzy 2004/2008: así reproduce IDT (docs/MOTOR_CALCULO.md §2.1)
    return mt.Tm_NN(s, nn_table=mt.DNA_NN3, dnac1=c["oligo_uM"] * 1000 / 4, dnac2=0,
                    Na=c["na_mM"], Mg=c["mg_mM"], dNTPs=c["dntp_mM"], saltcorr=7)


def _nn(dinuc):
    t = mt.DNA_NN3
    k = f"{dinuc}/{dinuc.translate(COMP)}"
    return t.get(k) or t[k[::-1]]


def dh_ds(s):
    dh, ds = mt.DNA_NN3["init"]
    for b in (s[0], s[-1]):
        h, e = mt.DNA_NN3["init_A/T" if b in "AT" else "init_G/C"]
        dh, ds = dh + h, ds + e
    for i in range(len(s) - 1):
        h, e = _nn(s[i:i + 2])
        dh, ds = dh + h, ds + e
    return dh, ds


def estabilidad_3p(s):
    dh = ds = 0
    for i in range(len(s) - 5, len(s) - 1):
        h, e = _nn(s[i:i + 2])
        dh, ds = dh + h, ds + e
    return dh - 310.15 * ds / 1000


def mw(s):
    return sum(MW[b] for b in s) - 61.964  # - PO2H + H2 (sin fosfato 5')


def eps260(s):
    return sum(EPS_NN[s[i:i + 2]] for i in range(len(s) - 1)) - sum(EPS_1[b] for b in s[1:-1])


def gc_pct(s):
    return 100 * sum(b in "GC" for b in s) / len(s)


def poly_x(s):
    return max(len(m.group()) for m in re.finditer(r"(.)\1*", s))


def _p3(c):
    return dict(mv_conc=c["na_mM"], dv_conc=c["mg_mM"], dntp_conc=c["dntp_mM"], dna_conc=c["oligo_uM"] * 1000)


def _lineal(s):
    return {"hebras": [{"seq": s, "k": list(range(len(s)))}], "pares": []}


def _dimero_3d(lineas):
    """Las 4 líneas de thal → hebras con su columna k (posición en la hélice) y pares de bases.
    Arriba: oligo A 5'→3'. Abajo: oligo B escrito 3'→5', como lo dibuja IDT."""
    filas = [x.split("\t", 1)[1] for x in lineas]
    ancho = max(map(len, filas))
    filas = [f.ljust(ancho) for f in filas]
    arriba, abajo, pares = {"seq": "", "k": []}, {"seq": "", "k": []}, []
    for k in range(ancho):
        a = next((f[k] for f in filas[:2] if f[k] not in " -"), None)
        b = next((f[k] for f in filas[2:] if f[k] not in " -"), None)
        for hebra, base in ((arriba, a), (abajo, b)):
            if base:
                hebra["seq"] += base
                hebra["k"].append(k)
        if filas[1][k] not in " -" and filas[2][k] not in " -":
            pares.append([len(arriba["seq"]) - 1, len(abajo["seq"]) - 1])
    return {"hebras": [arriba, abajo], "pares": pares}


def _horquilla_3d(s, marcas):
    """'/' = brazo 5' apareado, '\\' = brazo 3'. El brazo 3' se dobla bajo el 5' (escrito 3'→5')."""
    abre = [i for i, ch in enumerate(marcas) if ch == "/"]
    cierra = [i for i, ch in enumerate(marcas) if ch == "\\"]
    if not abre or len(abre) != len(cierra):
        return _lineal(s)
    pares_seq = list(zip(abre, reversed(cierra)))
    suma = abre[0] + cierra[-1]
    m = (abre[-1] + cierra[0]) // 2 + 1  # mitad del lazo
    js = list(range(len(s) - 1, m - 1, -1))
    pos = {j: t for t, j in enumerate(js)}
    return {"hebras": [{"seq": s[:m], "k": list(range(m))},
                       {"seq": "".join(s[j] for j in js), "k": [suma - j for j in js]}],
            "pares": [[i, pos[j]] for i, j in pares_seq], "lazo": True}


def horquilla(s, c, temp_c=25):
    r = primer3.calc_hairpin(s, temp_c=temp_c, output_structure=True, **_p3(c))
    if not r.structure_found:
        return {"dg": 0.0, "tm": None, "ascii": [], **_lineal(s)}
    lineas = r.ascii_structure_lines
    return {"dg": round(r.dg / 1000, 2), "tm": round(r.tm, 1), "ascii": lineas,
            **_horquilla_3d(s, lineas[0].split("\t", 1)[1])}


def dimero(a, b, c, temp_c=37):
    kw = dict(temp_c=temp_c, output_structure=True, **_p3(c))
    r = primer3.calc_homodimer(a, **kw) if a == b else primer3.calc_heterodimer(a, b, **kw)
    if not r.structure_found:
        return {"dg": 0.0, "ascii": [], "hebras": [{"seq": a, "k": list(range(len(a)))}], "pares": []}
    return {"dg": round(r.dg / 1000, 2), "ascii": r.ascii_structure_lines, **_dimero_3d(r.ascii_structure_lines)}


def anclaje_3p(a, b, c, temp_c=37):
    return round(min(primer3.calc_end_stability(x, y, temp_c=temp_c, **_p3(c)).dg for x, y in ((a, b), (b, a))) / 1000, 2)


def evaluar(valores, criterios):
    niveles, banderas = {}, []
    for etiqueta, campo, verde, amarillo in criterios:
        v = valores.get(campo)
        if v is None:
            continue
        n = "verde" if verde[0] <= v <= verde[1] else "amarillo" if amarillo[0] <= v <= amarillo[1] else "rojo"
        niveles[campo] = n
        if n != "verde":
            banderas.append({"regla": etiqueta, "campo": campo, "valor": v, "nivel": n})
    return niveles, banderas


def informe(nombre, seq, cond):
    s = limpiar(seq)
    fab, pcr = cond["fabrica"], cond["pcr"]
    dh, ds = dh_ds(s)
    m, e = mw(s), eps260(s)
    hp, sd = horquilla(s, fab), dimero(s, s, pcr)
    o = {
        "nombre": nombre, "seq": s, "largo": len(s), "gc_pct": round(gc_pct(s), 1),
        "tm_fabrica": round(tm(s, fab), 1), "tm_pcr": round(tm(s, pcr), 1),
        "pm": round(m, 1), "eps260": round(e), "nmol_od": round(1e6 / e, 2), "ug_od": round(1e3 * m / e, 2),
        "dh": round(dh, 1), "ds": round(ds, 1), "dg_duplex_37": round(dh - 310.15 * ds / 1000, 1),
        "horquilla_dg": hp["dg"], "horquilla_tm": hp["tm"],
        "autodimero_dg_fabrica": dimero(s, s, fab)["dg"], "autodimero_dg": sd["dg"],
        "pinza_gc": sum(b in "GC" for b in s[-5:]), "poly_x": poly_x(s), "estabilidad_3p": round(estabilidad_3p(s), 2),
        "estructuras": {"horquilla": hp, "autodimero": sd},
    }
    o["niveles"], o["banderas"] = evaluar(o, CRITERIOS_OLIGO)
    return o


def matriz_dimeros(oligos, c):
    """Todos contra todos (diagonal = autodímero). Ordenado del peor al mejor."""
    seqs = [limpiar(o["seq"]) for o in oligos]
    n = len(seqs)
    matriz = [[0.0] * n for _ in range(n)]
    pares = []
    for i in range(n):
        for j in range(i, n):
            d = dimero(seqs[i], seqs[j], c)
            anc = anclaje_3p(seqs[i], seqs[j], c)
            matriz[i][j] = matriz[j][i] = d["dg"]
            campo = "autodimero_dg" if i == j else "hetero_dg"
            niv, _ = evaluar({campo: d["dg"], "anclaje_3p": anc}, CRITERIOS_OLIGO + CRITERIOS_PAR)
            nivel = "rojo" if "rojo" in niv.values() else "amarillo" if "amarillo" in niv.values() else "verde"
            pares.append({"a": oligos[i]["nombre"], "b": oligos[j]["nombre"], "dg": d["dg"],
                          "anclaje_3p": anc, "nivel": nivel, "estructura": d})
    pares.sort(key=lambda p: p["dg"])
    return {"nombres": [o["nombre"] for o in oligos], "matriz": matriz, "pares": pares}
