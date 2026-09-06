#!/usr/bin/env python3
import os
import glob
import csv

SEGMENTOS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]

LONGITUDES_ESPERADAS = {
    "PB2": (2200, 2400),
    "PB1": (2200, 2400),
    "PA":  (2100, 2300),
    "HA":  (1600, 1800),
    "NP":  (1400, 1600),
    "NA":  (1300, 1550),
    "M":   (900, 1100),
    "NS":  (750, 950),
}

def leer_fasta(path):
    registros = []
    header = None
    seq = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    registros.append((header, "".join(seq).upper()))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())

    if header is not None:
        registros.append((header, "".join(seq).upper()))

    return registros

def detectar_segmento(header):
    h = header.upper()

    if "_HA_" in h or h.endswith("_HA") or h == "HA":
        return "HA"
    elif "_NA_" in h or h.endswith("_NA") or h == "NA":
        return "NA"
    elif "_PB2" in h or h.endswith("PB2") or h == "PB2":
        return "PB2"
    elif "_PB1" in h or h.endswith("PB1") or h == "PB1":
        return "PB1"
    elif "_PA" in h or h.endswith("PA") or h == "PA":
        return "PA"
    elif "_NP" in h or h.endswith("NP") or h == "NP":
        return "NP"
    elif "_NS" in h or h.endswith("NS") or h == "NS":
        return "NS"
    elif "_MP" in h or h.endswith("MP") or h == "MP" or h == "M":
        return "M"
    else:
        return "NO_DETECTADO"

def detectar_subtipo(header):
    h = header.upper()

    if "H1" in h or "N1" in h:
        return "H1N1"
    elif "H3" in h or "N2" in h:
        return "H3N2"
    else:
        return "NO_DETECTADO"

def calidad_por_n(pct_n):
    if pct_n < 1:
        return "Buena"
    elif pct_n < 5:
        return "Aceptable"
    elif pct_n < 10:
        return "Dudosa"
    else:
        return "Mala"

def estado_segmento(segmento, longitud, pct_n):
    if segmento not in LONGITUDES_ESPERADAS:
        return "Revisar"

    min_len, max_len = LONGITUDES_ESPERADAS[segmento]

    if longitud < min_len:
        return "Corto/incompleto"
    elif longitud > max_len:
        return "Longitud inesperada"
    elif pct_n >= 10:
        return "Muchos N"
    elif pct_n >= 5:
        return "Revisar por N"
    else:
        return "Correcto"

def main():
    fastas = sorted(
        glob.glob("*.fasta") +
        glob.glob("*.fa") +
        glob.glob("*.fas")
    )

    if not fastas:
        print("No se encontraron archivos FASTA en esta carpeta.")
        return

    filas_segmentos = []
    resumen = []

    for fasta in fastas:
        muestra = os.path.basename(fasta)
        muestra = muestra.replace(".consensus.fasta", "")
        muestra = muestra.replace(".fasta", "")
        muestra = muestra.replace(".fa", "")
        muestra = muestra.replace(".fas", "")

        registros = leer_fasta(fasta)

        segmentos_presentes = {}
        subtipo_muestra = "NO_DETECTADO"
        total_n = 0
        total_longitud = 0

        for header, seq in registros:
            segmento = detectar_segmento(header)
            subtipo = detectar_subtipo(header)

            if subtipo != "NO_DETECTADO":
                subtipo_muestra = subtipo

            longitud = len(seq)
            n_count = seq.count("N")
            pct_n = round((n_count / longitud) * 100, 2) if longitud > 0 else 0

            total_n += n_count
            total_longitud += longitud

            if segmento != "NO_DETECTADO":
                segmentos_presentes[segmento] = {
                    "longitud": longitud,
                    "n": n_count,
                    "pct_n": pct_n,
                    "estado": estado_segmento(segmento, longitud, pct_n)
                }

            filas_segmentos.append({
                "muestra": muestra,
                "archivo": fasta,
                "subtipo": subtipo if subtipo != "NO_DETECTADO" else subtipo_muestra,
                "segmento_detectado": segmento,
                "header_original": header,
                "longitud_nt": longitud,
                "N": n_count,
                "porcentaje_N": pct_n,
                "calidad_N": calidad_por_n(pct_n),
                "estado": estado_segmento(segmento, longitud, pct_n)
            })

        fila_resumen = {
            "muestra": muestra,
            "archivo": fasta,
            "subtipo": subtipo_muestra,
            "numero_segmentos_detectados": len([s for s in SEGMENTOS if s in segmentos_presentes]),
            "segmentos_presentes": ",".join([s for s in SEGMENTOS if s in segmentos_presentes]),
            "segmentos_ausentes": ",".join([s for s in SEGMENTOS if s not in segmentos_presentes]),
            "longitud_total": total_longitud,
            "N_total": total_n,
            "porcentaje_N_total": round((total_n / total_longitud) * 100, 2) if total_longitud > 0 else 0,
        }

        for seg in SEGMENTOS:
            if seg in segmentos_presentes:
                fila_resumen[f"{seg}_presente"] = "Sí"
                fila_resumen[f"{seg}_longitud"] = segmentos_presentes[seg]["longitud"]
                fila_resumen[f"{seg}_N"] = segmentos_presentes[seg]["n"]
                fila_resumen[f"{seg}_porcentaje_N"] = segmentos_presentes[seg]["pct_n"]
                fila_resumen[f"{seg}_estado"] = segmentos_presentes[seg]["estado"]
            else:
                fila_resumen[f"{seg}_presente"] = "No"
                fila_resumen[f"{seg}_longitud"] = ""
                fila_resumen[f"{seg}_N"] = ""
                fila_resumen[f"{seg}_porcentaje_N"] = ""
                fila_resumen[f"{seg}_estado"] = "Ausente"

        ha_ok = fila_resumen["HA_presente"] == "Sí" and fila_resumen["HA_estado"] == "Correcto"
        na_ok = fila_resumen["NA_presente"] == "Sí" and fila_resumen["NA_estado"] == "Correcto"

        if ha_ok and na_ok:
            fila_resumen["estado_general"] = "Apta para análisis HA/NA"
        elif fila_resumen["HA_presente"] == "Sí" or fila_resumen["NA_presente"] == "Sí":
            fila_resumen["estado_general"] = "Revisar HA/NA"
        else:
            fila_resumen["estado_general"] = "No apta para HA/NA"

        resumen.append(fila_resumen)

    with open("tabla_calidad_por_segmento.tsv", "w", newline="", encoding="utf-8") as out:
        campos = [
            "muestra",
            "archivo",
            "subtipo",
            "segmento_detectado",
            "header_original",
            "longitud_nt",
            "N",
            "porcentaje_N",
            "calidad_N",
            "estado"
        ]
        writer = csv.DictWriter(out, fieldnames=campos, delimiter="\t")
        writer.writeheader()
        writer.writerows(filas_segmentos)

    campos_resumen = list(resumen[0].keys())

    with open("tabla_resumen_por_muestra.tsv", "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=campos_resumen, delimiter="\t")
        writer.writeheader()
        writer.writerows(resumen)

    print("Análisis terminado.")
    print("Archivos generados:")
    print(" - tabla_calidad_por_segmento.tsv")
    print(" - tabla_resumen_por_muestra.tsv")

if __name__ == "__main__":
    main()
