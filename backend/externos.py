"""Conexiones externas: NCBI (bajar secuencias y BLAST de especificidad) e IDT (OligoAnalyzer oficial).

Variables de entorno (en Vercel: Settings → Environment Variables; nunca en el código):
  NCBI_EMAIL, NCBI_API_KEY                              opcionales; NCBI las pide para identificar la herramienta
  IDT_CLIENT_ID, IDT_CLIENT_SECRET, IDT_USUARIO, IDT_CONTRASENA   cuenta IDT con acceso a la API SciTools Plus
  PRIMERFLOW_CLAVE                                      clave que pide la página para usar IDT (que nadie más gaste tus consultas)
"""
import base64
import hmac
import json
import os
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
IDT = "https://www.idtdna.com"


def _pedir(url, datos=None, cabeceras=None, timeout=45):
    try:
        req = urllib.request.Request(url, data=datos, headers=cabeceras or {})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode()
    except urllib.error.HTTPError as e:
        raise ValueError(f"{url.split('/')[2]} respondió {e.code}: {e.read().decode(errors='replace')[:200]}") from None
    except (urllib.error.URLError, TimeoutError) as e:
        raise ValueError(f"No se pudo conectar con {url.split('/')[2]}: {e}") from None


def _ncbi():
    p = {"tool": "primerflow"}
    for clave, var in (("email", "NCBI_EMAIL"), ("api_key", "NCBI_API_KEY")):
        if os.environ.get(var):
            p[clave] = os.environ[var]
    return p


# ---------- NCBI: secuencias ----------
# «CP003967.2», «CP003967.2:1279779-1281089» o «CP003967.2:c1281089-1279779» (hebra complementaria, como en NCBI)
ACCESION = re.compile(r"^([A-Z]+_?[A-Z]*\d+(?:\.\d+)?)(?::(c?)(\d+)-(\d+))?$", re.I)


def ncbi_secuencias(consulta, maximo=30):
    """Accesiones (una por línea, con rango opcional) o una búsqueda libre de nuccore → texto FASTA."""
    lineas = [x.strip() for x in consulta.splitlines() if x.strip()]
    if not lineas:
        raise ValueError("Escribe una búsqueda o accesiones de NCBI")
    if all(ACCESION.match(x) for x in lineas):
        partes = []
        for x in lineas[:40]:  # ponytail: una consulta por rango a ≤3/s (límite de NCBI sin api key); ~40 caben en 60 s
            acc, comp, a, b = ACCESION.match(x).groups()
            p = {"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text"}
            if a:
                p |= {"seq_start": min(int(a), int(b)), "seq_stop": max(int(a), int(b)), "strand": 2 if comp else 1}
            partes.append(_pedir(EUTILS + "efetch.fcgi?" + urlencode(p | _ncbi())))
            time.sleep(0.34)
        return "".join(partes)
    termino = " ".join(lineas)
    if "[slen]" not in termino.lower():
        termino += " AND 100:20000[SLEN]"  # sin genomas completos
    termino += ' NOT "wgs master"[Properties]'  # los maestros de proyectos WGS no traen secuencia
    r = json.loads(_pedir(EUTILS + "esearch.fcgi?" + urlencode(
        {"db": "nuccore", "term": termino, "retmax": maximo, "retmode": "json"} | _ncbi())))
    ids = r["esearchresult"]["idlist"]
    if not ids:
        raise ValueError(f"NCBI no encontró nada para: {termino}")
    return _pedir(EUTILS + "efetch.fcgi", urlencode(
        {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"} | _ncbi()).encode())


# ---------- NCBI: BLAST de especificidad ----------
def blast_enviar(oligos, excluir=()):
    """Manda los primers a BLAST (blastn-short contra core_nt), sin las especies del panel → RID para preguntar después."""
    p = {"CMD": "Put", "PROGRAM": "blastn", "DATABASE": "core_nt", "MEGABLAST": "off", "WORD_SIZE": 7,
         "EXPECT": 1000, "HITLIST_SIZE": 100, "NUCL_REWARD": 1, "NUCL_PENALTY": -3, "GAPCOSTS": "5 2", "FILTER": "F",
         "QUERY": "".join(f">{o['nombre']}\n{o['seq']}\n" for o in oligos)}
    if excluir:  # «all[filter]» da el universo del que se resta (ENTREZ_QUERY no admite empezar con NOT)
        p["ENTREZ_QUERY"] = "all[filter] NOT (" + " OR ".join(f'"{e}"[Organism]' for e in excluir) + ")"
    html = _pedir(BLAST, urlencode(p | _ncbi()).encode())
    rid, rtoe = re.search(r"RID = (\w+)", html), re.search(r"RTOE = (\d+)", html)
    if not rid:
        raise ValueError("BLAST no aceptó la búsqueda")
    return {"rid": rid[1], "espera_s": int(rtoe[1]) if rtoe else 60}


def _sitios(reporte, max_mm):
    """Lugares donde el primer pega con su extremo 3′ completo y ≤ max_mm desajustes (huecos y 5′ sin alinear cuentan)."""
    b = reporte["report"]["results"]["search"]
    nombre, largo = b["query_title"], b["query_len"]
    for h in b.get("hits", []):
        d = h["description"][0]
        for s in h["hsps"]:
            mm = (s["query_from"] - 1) + (s["align_len"] - s["identity"])
            if s["query_to"] != largo or mm > max_mm:
                continue
            mas = s.get("hit_strand", "Plus") == "Plus"
            cinco = s["hit_from"] - (s["query_from"] - 1) if mas else s["hit_from"] + (s["query_from"] - 1)
            yield {"primer": nombre, "accesion": d["accession"], "especie": d.get("sciname", ""), "titulo": d.get("title", ""),
                   "hebra": "+" if mas else "-", "cinco": cinco, "tres": s["hit_to"], "desajustes": mm}


def blast_resultado(rid, pares, max_mm=3, tamano_max=4000):
    info = _pedir(BLAST + "?" + urlencode({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid} | _ncbi()))
    estado = (re.search(r"Status=(\w+)", info) or [None, "UNKNOWN"])[1]
    if estado == "WAITING":
        return {"estado": "esperando"}
    if estado != "READY":
        raise ValueError(f"BLAST {rid}: {estado}. La búsqueda falló o expiró (NCBI la guarda ~24 h); vuelve a enviarla.")
    if "ThereAreHits=yes" not in info:
        return {"estado": "listo", "sitios": [], "productos": []}
    datos = json.loads(_pedir(BLAST + "?" + urlencode({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2_S"} | _ncbi())))
    sitios = [s for r in datos["BlastOutput2"] for s in _sitios(r, max_mm)]
    return {"estado": "listo", "sitios": sitios, "productos": _productos(sitios, pares, tamano_max)}


def _productos(sitios, pares, tamano_max):
    # Producto: un sitio en + y otro en − del mismo registro, mirándose, a ≤ tamano_max (F-R, R-F, F-F o R-R del par).
    productos = []
    for par in pares:
        mios = [s for s in sitios if s["primer"] in (par["F"], par["R"])]
        for a in (s for s in mios if s["hebra"] == "+"):
            for b in (s for s in mios if s["hebra"] == "-" and s["accesion"] == a["accesion"]):
                tam = b["cinco"] - a["cinco"] + 1
                if a["tres"] < b["tres"] and 0 < tam <= tamano_max:
                    productos.append({"par": par["nombre"], "primers": f"{a['primer']} → {b['primer']}", "tamano": tam,
                                      "accesion": a["accesion"], "especie": a["especie"], "titulo": a["titulo"],
                                      "desajustes": a["desajustes"] + b["desajustes"]})
    return sorted(productos, key=lambda x: (x["desajustes"], x["tamano"]))


# ---------- IDT: OligoAnalyzer oficial (API SciTools Plus) ----------
def revisar_clave(clave):
    esperada = os.environ.get("PRIMERFLOW_CLAVE", "")
    if not esperada:
        raise ValueError("La conexión con IDT no está configurada en el servidor (falta PRIMERFLOW_CLAVE).")
    if not hmac.compare_digest(clave.encode(), esperada.encode()):
        raise PermissionError("Clave de acceso incorrecta.")


_token = {"valor": "", "vence": 0.0}


def _idt_token():
    if time.time() < _token["vence"]:
        return _token["valor"]
    v = {k: os.environ.get(k, "") for k in ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USUARIO", "IDT_CONTRASENA")}
    if not all(v.values()):
        raise ValueError("Faltan credenciales de IDT en el servidor: " + ", ".join(k for k, x in v.items() if not x))
    basic = base64.b64encode(f"{v['IDT_CLIENT_ID']}:{v['IDT_CLIENT_SECRET']}".encode()).decode()
    r = json.loads(_pedir(IDT + "/Identityserver/connect/token", urlencode(
        {"grant_type": "password", "username": v["IDT_USUARIO"], "password": v["IDT_CONTRASENA"], "scope": "test"}).encode(),
        {"Authorization": "Basic " + basic, "Content-Type": "application/x-www-form-urlencoded"}))
    _token.update(valor=r["access_token"], vence=time.time() + r.get("expires_in", 3600) - 60)
    return _token["valor"]


def _idt(ruta, cuerpo=None, consulta=None):
    url = IDT + "/restapi/v1/OligoAnalyzer/" + ruta + ("?" + urlencode(consulta) if consulta else "")
    return json.loads(_pedir(url, json.dumps(cuerpo).encode() if cuerpo else b"",
                             {"Authorization": "Bearer " + _idt_token(), "Content-Type": "application/json", "Accept": "application/json"}))


def _campo(d, *nombres):
    # ponytail: IDT no documenta su respuesta; se aceptan variantes de nombre y se devuelven los campos que llegaron
    k = {x.lower(): v for x, v in (d or {}).items()}
    return next((k[n.lower()] for n in nombres if n.lower() in k), None)


def _idt_oligo(o, cond):
    analizar = lambda c: _idt("Analyze", {"Sequence": o["seq"], "NaConc": c["na_mM"], "MgConc": c["mg_mM"],
                                          "dNTPsConc": c["dntp_mM"], "OligoConc": c["oligo_uM"], "NucleotideType": "DNA"})
    f = cond["fabrica"]
    fab, pcr = analizar(f), analizar(cond["pcr"])
    horquillas = _idt("Hairpin", {"Sequence": o["seq"], "NaConc": f["na_mM"], "MgConc": f["mg_mM"], "FoldingTemp": 25, "NucleotideType": "DNA"})
    dimeros = _idt("SelfDimer", consulta={"primary": o["seq"]})
    h = min(horquillas or [{}], key=lambda x: _campo(x, "deltaG") or 0)
    return {
        "nombre": o["nombre"],
        "tm_fabrica": _campo(fab, "MeltTemp"), "tm_pcr": _campo(pcr, "MeltTemp"),
        "pm": _campo(fab, "MolecularWeight", "MolWeight"), "eps260": _campo(fab, "ExtinctionCoefficient", "ExtCoefficient"),
        "nmol_od": _campo(fab, "NmoleOD", "NmolePerOD"), "ug_od": _campo(fab, "MicrogramOD", "UgOD", "MicrogramPerOD", "ugPerOD"),
        "horquilla_dg": _campo(h, "deltaG"), "horquilla_tm": _campo(h, "thermo", "Tm", "MeltTemp"),
        "autodimero_dg_fabrica": min((_campo(x, "DeltaG") for x in dimeros or []), default=None),
        "campos_idt": sorted(fab),
    }


def idt_oligos(oligos, cond):
    _idt_token()  # una sola vez antes de repartir el trabajo
    with ThreadPoolExecutor(max_workers=4) as ex:
        return list(ex.map(lambda o: _idt_oligo(o, cond), oligos))


if __name__ == "__main__":  # comprobación sin red: sitios y productos de un reporte BLAST mínimo
    rep = lambda q, hs: {"report": {"results": {"search": {"query_title": q, "query_len": 20, "hits": [
        {"description": [{"accession": "X1", "sciname": "E. coli", "title": "t"}], "hsps": hs}]}}}}
    F = rep("F", [{"query_from": 1, "query_to": 20, "align_len": 20, "identity": 19, "hit_from": 100, "hit_to": 119, "hit_strand": "Plus"}])
    R = rep("R", [{"query_from": 2, "query_to": 20, "align_len": 19, "identity": 19, "hit_from": 330, "hit_to": 312, "hit_strand": "Minus"},
                  {"query_from": 1, "query_to": 18, "align_len": 18, "identity": 18, "hit_from": 900, "hit_to": 883, "hit_strand": "Minus"}])
    s = [x for r in (F, R) for x in _sitios(r, 3)]
    assert len(s) == 2 and s[1]["cinco"] == 331 and s[1]["desajustes"] == 1, s  # el 2.º sitio de R no cubre el 3′
    p = _productos(s, [{"nombre": "par", "F": "F", "R": "R"}], 4000)
    assert [(x["primers"], x["tamano"]) for x in p] == [("F → R", 232)], p
    assert _productos(s, [{"nombre": "par", "F": "F", "R": "R"}], 200) == []
    assert ACCESION.match("CP003967.2:c1281089-1279779").groups() == ("CP003967.2", "c", "1281089", "1279779")
    assert ACCESION.match("NZ_CP012345.1") and not ACCESION.match("Acinetobacter baumannii basC")
    print("ok")
