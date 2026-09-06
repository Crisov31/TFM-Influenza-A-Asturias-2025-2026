#!/usr/bin/env python3
import os
import glob
import csv

SEGMENTOS_INTERES = ["HA", "NA"]

LONGITUDES_ESPERADAS = {
    "HA": (1600, 1800),
    "NA": (1300, 1550),
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
                seq.append(line)

    if header is not None:
        registros.append((header, "".join(seq).upper()))

    return registros

def detectar_segmento(header):
    h = header.upper()

    if "_HA_" in h or h.startswith("A_HA") or h == "HA":
        return "HA"
    elif "_NA_" in h or h.startswith("A_NA") or h == "NA":
        return "NA"
    elif "_PB2" in h or h.endswith("PB2"):
        return "PB2"
    elif "_PB1" in h or h.endswith("PB1"):
        return "PB1"
    elif "_PA" in h or h.endswith("PA"):
        return "PA"
    elif "_NP" in h or h.endswith("NP"):
        return "NP"
    elif "_NS" in h or h.endswith("NS"):
        return "NS"
    elif "_MP" in h or h.endswith("MP"):
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

def limpiar_nombre_muestra(nombre_archivo):
    nombre = os.path.basename(nombre_archivo)
    for ext in [".consensus.fasta", ".fasta", ".fa", ".fas"]:
        nombre = nombre.replace(ext, "")
    return nombre

def posiciones_bases_reales(seq):
    """
    Devuelve primera y última posición con base real A/C/G/T.
    Posiciones en formato 1-based.
    """
    bases_reales = {"A", "C", "G", "T"}

    primera = None
    ultima = None

    for i, base in enumerate(seq):
        if base in bases_reales:
            primera = i
            break

    for i in range(len(seq) - 1, -1, -1):
        if seq[i] in bases_reales:
            ultima = i
            break

    if primera is None or ultima is None:
        return None, None

    return primera + 1, ultima + 1

def clasificar_calidad(cobertura_pct, pct_n_internas):
    if cobertura_pct >= 90 and pct_n_internas < 1:
        return "Buena"
    elif cobertura_pct >= 80 and pct_n_internas < 5:
        return "Aceptable"
    elif cobertura_pct >= 60 and pct_n_internas < 20:
        return "Dudosa"
    else:
        return "Excluir"

def main():
    fastas = sorted(
        glob.glob("*.fasta") +
        glob.glob("*.fa") +
        glob.glob("*.fas")
    )

    if not fastas:
        print("No se encontraron archivos FASTA.")
        return

    os.makedirs("segmentos_filtrados", exist_ok=True)

    salidas = {
        "HA_H1N1": open("segmentos_filtrados/HA_H1N1.fasta", "w", encoding="utf-8"),
        "HA_H3N2": open("segmentos_filtrados/HA_H3N2.fasta", "w", encoding="utf-8"),
        "NA_H1N1": open("segmentos_filtrados/NA_H1N1.fasta", "w", encoding="utf-8"),
        "NA_H3N2": open("segmentos_filtrados/NA_H3N2.fasta", "w", encoding="utf-8"),
    }

    filas = []

    for fasta in fastas:
        muestra = limpiar_nombre_muestra(fasta)
        registros = leer_fasta(fasta)

        subtipo_muestra = "NO_DETECTADO"
        for header, seq in registros:
            subtipo_tmp = detectar_subtipo(header)
            if subtipo_tmp != "NO_DETECTADO":
                subtipo_muestra = subtipo_tmp
                break

        for header, seq in registros:
            segmento = detectar_segmento(header)

            if segmento not in SEGMENTOS_INTERES:
                continue

            subtipo = detectar_subtipo(header)
            if subtipo == "NO_DETECTADO":
                subtipo = subtipo_muestra

            longitud_total = len(seq)
            n_total = seq.count("N")

            primera, ultima = posiciones_bases_reales(seq)

            if primera is None:
                longitud_cubierta = 0
                region_cubierta = ""
                n_internas = 0
                pct_n_internas = 100
                cobertura_pct = 0
            else:
                region_cubierta = seq[primera-1:ultima]
                longitud_cubierta = len(region_cubierta)
                n_internas = region_cubierta.count("N")
                pct_n_internas = round((n_internas / longitud_cubierta) * 100, 2) if longitud_cubierta > 0 else 100

                min_esperada, max_esperada = LONGITUDES_ESPERADAS[segmento]
                cobertura_pct = round((longitud_cubierta / max_esperada) * 100, 2)

            calidad = clasificar_calidad(cobertura_pct, pct_n_internas)

            utilizable = "Sí" if calidad in ["Buena", "Aceptable", "Dudosa"] else "No"

            filas.append({
                "muestra": muestra,
                "archivo": fasta,
                "subtipo": subtipo,
                "segmento": segmento,
                "header_original": header,
                "longitud_total": longitud_total,
                "N_total": n_total,
                "primera_base_real": primera if primera is not None else "",
                "ultima_base_real": ultima if ultima is not None else "",
                "longitud_cubierta": longitud_cubierta,
                "cobertura_pct_aprox": cobertura_pct,
                "N_internas": n_internas,
                "porcentaje_N_internas": pct_n_internas,
                "calidad": calidad,
                "utilizable": utilizable
            })

            if utilizable == "Sí" and subtipo in ["H1N1", "H3N2"]:
                clave = f"{segmento}_{subtipo}"

                if clave in salidas:
                    nuevo_header = f">{muestra}|{segmento}|{subtipo}|calidad={calidad}|cobertura={cobertura_pct}|Nint={pct_n_internas}"
                    salidas[clave].write(nuevo_header + "\n")
                    salidas[clave].write(seq + "\n")

    for f in salidas.values():
        f.close()

    with open("tabla_calidad_cobertura_real.tsv", "w", newline="", encoding="utf-8") as out:
        campos = [
            "muestra",
            "archivo",
            "subtipo",
            "segmento",
            "header_original",
            "longitud_total",
            "N_total",
            "primera_base_real",
            "ultima_base_real",
            "longitud_cubierta",
            "cobertura_pct_aprox",
            "N_internas",
            "porcentaje_N_internas",
            "calidad",
            "utilizable"
        ]

        writer = csv.DictWriter(out, fieldnames=campos, delimiter="\t")
        writer.writeheader()
        writer.writerows(filas)

    print("Análisis terminado.")
    print("Archivos generados:")
    print(" - tabla_calidad_cobertura_real.tsv")
    print(" - segmentos_filtrados/HA_H1N1.fasta")
    print(" - segmentos_filtrados/HA_H3N2.fasta")
    print(" - segmentos_filtrados/NA_H1N1.fasta")
    print(" - segmentos_filtrados/NA_H3N2.fasta")

if __name__ == "__main__":
    main()
