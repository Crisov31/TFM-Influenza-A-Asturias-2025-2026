#!/usr/bin/env python3
import csv

# Posiciones clave (numeración NA estándar)
MUTACIONES = {
    "H1N1": {
        119: ["V", "G", "D", "A"],
        136: ["K"],
        223: ["R", "K", "V"],
        275: ["Y"],
        295: ["S"],
    },
    "H3N2": {
        119: ["V", "I"],
        136: ["K"],
        222: ["V", "T", "R"],
        292: ["K"],
        294: ["S"],
    }
}

def leer_fasta(path):
    seqs = {}
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.upper())
        if header:
            seqs[header] = "".join(seq)

    return seqs

def encontrar_referencia(seqs):
    for h in seqs:
        if "REF" in h.upper():
            return h
    raise Exception("No se encontró referencia (REF)")

def mapear_posiciones_ref(seq_ref):
    mapa = {}
    pos = 0

    for i, aa in enumerate(seq_ref):
        if aa != "-":
            pos += 1
            mapa[pos] = i

    return mapa

def analizar(alineado, subtipo, salida):
    seqs = leer_fasta(alineado)
    ref_name = encontrar_referencia(seqs)

    ref_seq = seqs[ref_name]
    mapa = mapear_posiciones_ref(ref_seq)

    filas = []

    for nombre, seq in seqs.items():
        if nombre == ref_name:
            continue

        estados = []

        for pos, aa_res in MUTACIONES[subtipo].items():
            if pos not in mapa:
                estado = "NO_CUBIERTO"
                aa = ""
            else:
                i = mapa[pos]
                aa = seq[i]

                if aa == "-":
                    estado = "NO_CUBIERTO"
                elif aa == "X":
                    estado = "NO_EVALUABLE"
                elif aa in aa_res:
                    estado = "MUTACION_RESISTENCIA"
                else:
                    estado = "SIN_MUTACION"

            estados.append(estado)

            filas.append({
                "muestra": nombre,
                "subtipo": subtipo,
                "posicion": pos,
                "aa_detectado": aa,
                "aa_resistencia": ",".join(aa_res),
                "estado": estado
            })

    with open(salida, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=filas[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(filas)

    print(f"Generado {salida}")

def main():
    analizar("NA_H1N1_aligned.fasta", "H1N1", "resistencias_NA_H1N1.tsv")
    analizar("NA_H3N2_aligned.fasta", "H3N2", "resistencias_NA_H3N2.tsv")

if __name__ == "__main__":
    main()
