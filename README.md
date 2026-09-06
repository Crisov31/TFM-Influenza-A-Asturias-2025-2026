# TFM – Influenza A Asturias 2025–2026

This repository contains the scripts and supplementary data generated for the Master's Thesis focused on the genomic surveillance of Influenza A viruses circulating in Asturias during the 2025–2026 season.

The analyses included sequence quality control, filtering of HA and NA segments, phylogenetic dataset preparation, and screening of neuraminidase (NA) amino acid substitutions associated with antiviral resistance.

## Repository structure

### `scripts/`

Python scripts used during sequence processing and analysis.

#### `control_calidad_fastas.py`

Performs an initial quality assessment of Influenza A consensus FASTA files.

The script:
- detects the eight genomic segments: PB2, PB1, PA, HA, NP, NA, M and NS;
- identifies the viral subtype when possible;
- calculates sequence length and percentage of ambiguous bases (N);
- evaluates the quality of each genomic segment;
- generates summary tables at segment and sample level.

Outputs:

- `tabla_calidad_por_segmento.tsv`
- `tabla_resumen_por_muestra.tsv`

#### `control_y_filtrado_HA_NA.py`

Performs quality control specifically for the HA and NA segments.

The script:
- identifies HA and NA sequences;
- calculates the covered region between the first and last real nucleotide;
- estimates sequence coverage;
- calculates the proportion of internal ambiguous bases;
- classifies sequences according to quality criteria;
- selects sequences considered suitable for subsequent analyses.

Sequences with at least 60% approximate coverage and less than 20% internal ambiguous bases can be retained for downstream analyses.

Outputs:

- `tabla_calidad_cobertura_real.tsv`
- filtered FASTA files for HA and NA of H1N1 and H3N2 viruses.

#### `filtrar_Ns.py`

Filters FASTA sequences according to the proportion of ambiguous nucleotides.

Sequences containing less than 20% ambiguous bases (N) are retained.

Usage:

```bash
python filtrar_Ns.py input.fasta output.fasta
```

#### `resistencias_desde_alineamiento.py`

Screens aligned neuraminidase protein sequences for amino acid substitutions associated with reduced susceptibility or resistance to neuraminidase inhibitors.

The script evaluates selected positions in H1N1 and H3N2 NA sequences using a reference sequence included in the alignment.

Outputs:

- `resistencias_NA_H1N1.tsv`
- `resistencias_NA_H3N2.tsv`

---

## Supplementary tables

The directory `supplementary_tables/` contains the tables generated during sequence quality control and antiviral resistance analyses.

### `tabla_calidad_por_segmento.tsv`

Quality assessment of each genomic segment detected in every sample.

### `tabla_resumen_por_muestra.tsv`

Summary of the genomic segments detected for each sample and their overall quality status.

### `tabla_calidad_cobertura_real.tsv`

Quality assessment of HA and NA sequences based on approximate coverage and internal ambiguous nucleotide content.

### `resistencias_NA_H1N1.tsv`

Results of the screening of H1N1 neuraminidase sequences for amino acid substitutions associated with antiviral resistance.

### `resistencias_NA_H3N2.tsv`

Results of the screening of H3N2 neuraminidase sequences for amino acid substitutions associated with antiviral resistance.

---

## Software

The analyses were performed using Python scripts together with standard bioinformatics tools used for sequence alignment, phylogenetic analysis and clade classification.

Main tools used during the study included:

- Python
- MAFFT
- IQ-TREE
- Nextclade
- iTOL
- FigTree

---

## Study context

The repository accompanies the Master's Thesis investigating the genomic characteristics of Influenza A viruses detected in Asturias during the 2025–2026 influenza season.

The genomic analyses focused primarily on the hemagglutinin (HA) and neuraminidase (NA) segments for phylogenetic characterization and antiviral resistance surveillance.

---

## Author

Cristina Ochoa Varela

Master
2025–2026
