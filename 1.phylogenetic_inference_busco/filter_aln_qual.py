import os
import subprocess
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict
import csv
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--pathout")
parser.add_argument("--fasta")
parser.add_argument("--stats")
args = parser.parse_args()

# The parameters used
input_dir = args.path
clean_dir = args.pathout
concat_file = args.fasta
report_file = args.stats

# The filters used
max_gap_pct = 5
min_aln_length = 2000
min_variable_sites = 5

os.makedirs(clean_dir, exist_ok=True)

def compute_alignment_stats(alignment):
    aln_len = alignment.get_alignment_length()
    n_seq = len(alignment)
    gap_count = sum(record.seq.count("-") for record in alignment)
    total_sites = aln_len * n_seq
    gap_percentage = (gap_count / total_sites) * 100 if total_sites > 0 else 0

    columns = [str(alignment[:, i]) for i in range(aln_len)]
    conserved = sum(1 for col in columns if len(set(col.replace("-", ""))) == 1)
    variable = sum(1 for col in columns if len(set(col.replace("-", ""))) > 1)

    return {
        "n_seq": n_seq,
        "aln_length": aln_len,
        "gap_percentage": round(gap_percentage, 2),
        "constant_sites": conserved,
        "variable_sites": variable
    }

# Step 1: Cleansing and stats
qualified_alignments = []
all_taxa = set()
alignment_dict = {}
c=0
with open(report_file, "w", newline="") as csvfile:
    fieldnames = ["gene", "n_seq", "aln_length", "gap_percentage", "constant_sites", "variable_sites"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for fname in os.listdir(input_dir):
        if not fname.endswith(".fasta"):
            continue
        c+=1
        print(c,fname)
        gene = fname.replace(".fasta", "")
        in_path = os.path.join(input_dir, fname)
        out_path = os.path.join(clean_dir, fname)

        subprocess.run(["trimal", "-in", in_path, "-out", out_path, "-automated1"], check=False)

        try:
            alignment = AlignIO.read(out_path, "fasta")
            stats = compute_alignment_stats(alignment)
            stats["gene"] = gene
            writer.writerow(stats)

            if (
                stats["gap_percentage"] <= max_gap_pct and
                stats["aln_length"] >= min_aln_length and
                stats["variable_sites"] >= min_variable_sites
            ):
                qualified_alignments.append((gene, alignment))
                all_taxa.update(rec.id for rec in alignment)
                alignment_dict[gene] = alignment

        except Exception as e:
            print(f"[!] Erreur dans {fname} : {e}")

# Step 2 : concatenation
if len(qualified_alignments) == 0:
    print("Aucun alignement ne passe les filtres. Aucune concaténation effectuée.")
else:
    all_taxa = sorted(all_taxa)
    supermatrix = defaultdict(str)

    for gene, alignment in qualified_alignments:
        gene_len = alignment.get_alignment_length()
        id_to_seq = {rec.id: str(rec.seq) for rec in alignment}

        for taxon in all_taxa:
            if taxon in id_to_seq:
                supermatrix[taxon] += id_to_seq[taxon]
            else:
                supermatrix[taxon] += "-" * gene_len

    # Writing the supermatrix
    with open(concat_file, "w") as out_f:
        for taxon in all_taxa:
            out_f.write(f">{taxon}\n{supermatrix[taxon]}\n")

    print(f"{len(qualified_alignments)} alignements retenus")
    print(f"Supermatrice écrite dans : {concat_file}")
print(f"Statistiques dans : {report_file}")

for gene, alignment in qualified_alignments:
    gene_len = alignment.get_alignment_length()

# Step 3 : Generating the partitions files 
partitions_file = concat_file.replace(".fasta", ".partitions.txt")

start = 1
with open(partitions_file, "w") as pf:
    for gene, alignment in qualified_alignments:
        gene_len = alignment.get_alignment_length()
        end = start + gene_len - 1
        pf.write(f"DNA, {gene} = {start}-{end}\n")
        start = end + 1

print(f"Partitions écrites dans : {partitions_file}")

