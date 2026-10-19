"""Verifica una tabla pegada (tus columnas, en cualquier orden) contra el motor, celda por celda.

Estados: ok = coincide · aprox = difiere por un modelo distinto y conocido (IDT con Mg²⁺, UNAFold) · no = no coincide."""
import re
import unicodedata

import motor

INF = float("inf")
NOTA_MG = ("IDT con Mg²⁺ se separa de las ecuaciones publicadas (Owczarzy 2008): −0.4 °C de media, hasta 0.7. "
           "Ninguna combinación de condiciones lo reproduce exacto; es una diferencia de implementación de IDT.")
NOTA_UNAFOLD = "IDT calcula horquillas con UNAFold (de pago); el motor usa thal de Primer3. Tolerancia de tu manual: ±1.5 kcal/mol."

# encabezado normalizado → (campo del informe, tolerancia "ok", tolerancia "aprox", nota si sólo es aprox, cómo se calcula)
COLUMNAS = {
    "largo": ("largo", 0, 0, "", "nº de bases"),
    "gc pct": ("gc_pct", 0.051, 0.051, "", "(G + C) / largo × 100"),
    "tm idt valores de fabrica": ("tm_fabrica", 0.11, 0.5, "", "NN SantaLucia 1998 + Owczarzy 2004, preset fábrica"),
    "tm idt condiciones pcr": ("tm_pcr", 0.11, 1.0, NOTA_MG, "NN SantaLucia 1998 + Owczarzy 2008, preset PCR"),
    "peso molecular g mol": ("pm", 0.11, 0.11, "", "Σ masas de IDT − 61.964"),
    "coef extincion l mol cm": ("eps260", 1, 1, "", "vecino más cercano (Cavaluzzi y Borer)"),
    "nmol por od260": ("nmol_od", 0.011, 0.011, "", "10⁶ / ε260"),
    "ug por od260": ("ug_od", 0.011, 0.011, "", "PM / ε260 × 1000"),
    "dh kcal mol entalpia": ("dh", 0.06, 0.06, "", "Σ ΔH vecino más cercano + iniciación"),
    "ds cal k mol entropia": ("ds", 0.06, 0.06, "", "Σ ΔS vecino más cercano + iniciación"),
    "horquilla dg fabrica 25c": ("horquilla_dg", 0.011, 1.5, NOTA_UNAFOLD, "thal (Primer3), 25 °C, preset fábrica"),
    "horquilla tm": ("horquilla_tm", 0.11, INF, NOTA_UNAFOLD + " La Tm de horquilla no es comparable entre motores.", "thal (Primer3)"),
    "autodimero dg fabrica": ("autodimero_dg_fabrica", 0.011, 0.011, "", "thal (Primer3), 37 °C, preset fábrica"),
    "autodimero dg pcr": ("autodimero_dg", 0.011, 0.011, "", "thal (Primer3), 37 °C, preset PCR"),
    "dg maximo duplex perfecto": ("dg_duplex_37", 0.06, 0.06, "", "ΔH − 310.15 × ΔS"),
}


def norm(t):
    t = unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", " ", t).strip()


def numero(v):
    try:
        return float(str(v).strip().replace("−", "-").replace(",", "."))
    except ValueError:
        return None


def celda(tuyo, calc, tol, tol_aprox, nota_aprox, como):
    t = numero(tuyo)
    if t is None:
        return {"tuyo": tuyo, "estado": ""}
    if calc is None:  # thal no encuentra horquilla estable
        estado = "aprox" if t >= 0 or tol_aprox == INF else "no"
        return {"tuyo": tuyo, "calculado": None, "estado": estado,
                "nota": "El motor no encuentra horquilla estable (ΔG ≥ 0). " + nota_aprox}
    d = round(calc - t, 2)
    if abs(d) <= tol:
        return {"tuyo": tuyo, "calculado": calc, "diferencia": d, "estado": "ok"}
    if calc == 0 and t >= 0 and "horquilla" in como + nota_aprox:
        return {"tuyo": tuyo, "calculado": calc, "diferencia": d, "estado": "aprox",
                "nota": "Equivalente: ninguno de los dos predice una horquilla estable (ΔG ≥ 0). " + nota_aprox}
    if abs(d) <= tol_aprox:
        return {"tuyo": tuyo, "calculado": calc, "diferencia": d, "estado": "aprox", "nota": nota_aprox}
    return {"tuyo": tuyo, "calculado": calc, "diferencia": d, "estado": "no",
            "nota": f"El valor correcto es {calc} ({como})." + (" " + nota_aprox if nota_aprox else "")}


def verificar_oligos(encabezados, filas, cond):
    claves = [norm(h) for h in encabezados]
    i_seq = next((i for i, k in enumerate(claves) if k.startswith("secuencia") or k in ("seq", "sequence")), None)
    if i_seq is None:
        raise ValueError("La tabla necesita una columna de secuencia (p. ej. «secuencia 5a3»).")
    salida, resumen = [], {"ok": 0, "aprox": 0, "no": 0}
    for fila in filas:
        fila = list(fila) + [""] * (len(encabezados) - len(fila))
        inf = motor.informe("x", fila[i_seq], cond)
        celdas = []
        for k, v in zip(claves, fila):
            if k in COLUMNAS:
                campo, tol, tol_a, nota, como = COLUMNAS[k]
                c = celda(v, inf[campo], tol, tol_a, nota, como)
                if c["estado"]:
                    resumen[c["estado"]] += 1
            else:
                c = {"tuyo": v, "estado": ""}
            celdas.append(c)
        salida.append(celdas)
    return {"encabezados": encabezados, "filas": salida, "resumen": resumen}
