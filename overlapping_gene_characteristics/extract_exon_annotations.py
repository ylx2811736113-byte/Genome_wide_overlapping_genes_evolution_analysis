import csv

# Define input and output file paths
transcript_file = (
    "D:/Work/overlapping_gene/results/animal_new/Pocillopora_verrucosa/protein_gene_info.txt_longest_tra.csv"
)

gff_file = (
    "D:/Work/overlapping_gene/GFF/animal_new/Pocillopora_verrucosa/Pocillopora_verrucosa_genomic.gff"
)

output_file = (
    "D:/Work/overlapping_gene/results/animal_new/Pocillopora_verrucosa/exon.txt"
)


def parse_attributes(attribute_string):
    """
    Parse the GFF attribute column into a dictionary.
    """
    attributes = {}

    for item in attribute_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value

    return attributes


# Build a transcript-to-chromosome lookup table
transcript_to_chr = {}

with open(
    transcript_file,
    "r",
    newline="",
    encoding="utf-8"
) as infile:

    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:

        transcript_to_chr[
            row["transcript_ID"]
        ] = row["transcript_chr_id"]


# Extract exon annotations
with open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as outfile:

    writer = csv.writer(outfile, delimiter="\t")

    writer.writerow([
        "chromosome",
        "start",
        "end",
        "score",
        "strand",
        "exon_ID",
        "transcript_ID",
    ])

    with open(
        gff_file,
        "r",
        encoding="utf-8"
    ) as infile:

        for line in infile:

            if line.startswith("#") or not line.strip():
                continue

            columns = line.rstrip().split("\t")

            if len(columns) != 9:
                continue

            (
                seqid,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = columns

            if feature != "exon":
                continue

            attr = parse_attributes(attributes)

            parent = attr.get("Parent", "")
            exon_id = attr.get("ID", "")

            transcript_id = (
                parent.split("rna-")[-1]
                if "rna-" in parent
                else parent
            )

            exon_id = (
                exon_id.split("exon-")[-1]
                if "exon-" in exon_id
                else exon_id
            )

            chromosome = transcript_to_chr.get(transcript_id)

            if chromosome is None:
                continue

            writer.writerow([
                chromosome,
                start,
                end,
                score,
                strand,
                exon_id,
                transcript_id,
            ])

print(f"Results have been saved to {output_file}")
