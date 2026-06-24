import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--pathout")
parser.add_argument("--pathoutdir")
args = parser.parse_args()

import os
import subprocess
from Bio import AlignIO
import csv

input_dir = args.path
output_dir = args.pathoutdir
report_file = args.pathout

os.makedirs(output_dir, exist_ok=True)

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

with open(report_file, "w", newline="") as csvfile:
    fieldnames = ["gene", "n_seq", "aln_length", "gap_percentage", "constant_sites", "variable_sites"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for fname in os.listdir(input_dir):
        if not fname.endswith(".fasta"):
            continue
        print(fname)
        gene = fname.replace(".fasta", "")
        in_path = os.path.join(input_dir, fname)
        out_path = os.path.join(output_dir, fname)

        # Clean with trimAl
        subprocess.run(["trimal", "-in", in_path, "-out", out_path, "-automated1"])

        try:
            alignment = AlignIO.read(out_path, "fasta")
            stats = compute_alignment_stats(alignment)
            stats["gene"] = gene
            writer.writerow(stats)
        except Exception as e:
            print(f"Erreur avec {fname} : {e}")
