# Overlapping Gene Identification and Classification

This directory contains the scripts used to identify and classify overlapping protein-coding genes from genome annotation (GFF/GFF3) files. Because the annotation structures differ substantially among animals, plants, and prokaryotes, separate workflows are provided for each group.

## Workflow

The pipeline consists of the following scripts:

| Script | Description |
|---------|-------------|
| `extract_gff_information.py` | Extract protein-coding gene, mRNA, and CDS annotations from complex animal GFF/GFF3 files. |
| `extract_longest_transcript.py` | Identify the longest transcript for each protein-coding gene. |
| `extract_cds.py` | Extract CDS annotations associated with the longest transcript. |
| `extract_cds_prokaryote.py` | Extract CDS annotations from prokaryotic GFF/GFF3 files. |
| `identify_gene_coordinates_via_longest_transcript_cds.py` | Determine representative gene coordinates using the CDS regions of the longest transcript. |
| `identify_overlapping_genes.py` | Identify overlapping genes and classify overlap types. |
| `main_script_animal.py` | Execute the complete workflow for animal genomes. |
| `main_script_plant.py` | Execute the complete workflow for plant genomes. |
| `main_script_prokaryote.py` | Execute the complete workflow for prokaryotic genomes. |

## Workflow for Different Taxonomic Groups

### Animal genomes

Animal GFF/GFF3 files typically contain complex hierarchical annotations (e.g., gene–mRNA–CDS relationships). Therefore, the complete workflow requires all of the following scripts:

- `extract_gff_information.py`
- `extract_longest_transcript.py`
- `extract_cds.py`
- `identify_gene_coordinates_via_longest_transcript_cds.py`
- `identify_overlapping_genes.py`
- `main_script_animal.py`

---

### Plant genomes

Plant GFF/GFF3 annotations are relatively simpler and do not require the additional preprocessing performed for animal annotations. The plant workflow consists of:

- `extract_longest_transcript.py`
- `extract_cds.py`
- `identify_gene_coordinates_via_longest_transcript_cds.py`
- `identify_overlapping_genes.py`
- `main_script_plant.py`

---

### Prokaryotic genomes

Prokaryotic genome annotations are generally simple and usually contain only CDS features without transcript isoforms. Consequently, only three scripts are required:

- `extract_cds_prokaryote.py`
- `identify_overlapping_genes.py`
- `main_script_prokaryote.py`

## Notes

- The three `main_script_*.py` programs serve as workflow entry points and automatically execute the required scripts for each taxonomic group.
- Representative gene coordinates are defined using the CDS coordinates of the longest transcript for eukaryotic genomes (animals and plants).
- Overlapping genes are subsequently identified and classified according to their genomic coordinates and transcriptional orientations.
