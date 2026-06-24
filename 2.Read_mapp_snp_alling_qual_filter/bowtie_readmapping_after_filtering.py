import os
import argparse
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Generate and submit SLURM jobs for Bowtie2 mapping.")
    parser.add_argument("--pathin", required=True, help="Base path/prefix for input fastq files")
    parser.add_argument("--ref", required=True, help="Directory containing reference fasta files")
    parser.add_argument("--script_dir", default="./scripts", help="Directory to save .sh scripts")
    parser.add_argument("--pathout", default="./bam", help="Directory for output BAM files")
    parser.add_argument("--log_dir", default="./logs", help="Directory for SLURM logs")
    parser.add_argument("--sensitive", default="sensitive", help="Bowtie2 sensitivity (e.g., very-sensitive)")
    args = parser.parse_args()

    # Convert to Path objects for easier manipulation
    path_in = Path(args.pathin)
    ref_dir = Path(args.ref)
    script_dir = Path(args.script_dir)
    path_out = Path(args.pathout)
    log_dir = Path(args.log_dir)

    # Ensure directories exist
    for d in [script_dir, path_out, log_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # Metadata for naming
    ind_name = path_in.name

    # Detect naming format
    if os.path.exists(f"{path_in}_FW_paired.fq.gz"):
        reads1 = f"{path_in}_FW_paired.fq.gz"
        reads2 = f"{path_in}_RV_paired.fq.gz"
    elif os.path.exists(f"{path_in}.R1.fastq.gz"):
        reads1 = f"{path_in}.R1.fastq.gz"
        reads2 = f"{path_in}.R2.fastq.gz"
    elif os.path.exists(f"{path_in}_R1_merged.fastq.gz"):
        reads1 = f"{path_in}_R1_merged.fastq.gz"
        reads2 = f"{path_in}_R2_merged.fastq.gz"    
    else:
        print(f"No reads found for {path_in}, skipping.")
        return

    # Loop over reference genomes
    for ref_file in ref_dir.glob("*"):
        if ref_file.suffix not in ['.fasta', '.fa']:
            continue

        # Clean reference name (removing common assembly suffixes)
        ref_name = ref_file.name.replace('_assembly.fasta', '').replace('-scaffolds.fa', '').replace(ref_file.suffix, '')

        # Define file paths
        job_tag = f"{ind_name}_{ref_name}_{args.sensitive}"
        print(job_tag)
        script_path = script_dir / f"{job_tag}.sh"
        print(script_path)
        bowtie_index = ref_file.with_suffix('')  # Assumes index prefix matches fasta prefix
        output_bam = path_out / f"{job_tag}.bam"

        # Create the SLURM script using a multi-line f-string
        slurm_content = f"""#!/bin/bash
#SBATCH -J {job_tag}
#SBATCH -o {log_dir}/{job_tag}.out
#SBATCH -e {log_dir}/{job_tag}.err
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=8
module load bowtie2/2.3.4.3
module load samtools/1.9
module load fastp/1.0.1
# Step 1 — Trim reads with fastp
fastp \\
    --in1 {reads1} \\
    --in2 {reads2} \\
    --out1 {path_out}/{ind_name}_trimmed_R1.fastq.gz \\
    --out2 {path_out}/{ind_name}_trimmed_R2.fastq.gz \\
    --thread 8 \\
    --qualified_quality_phred 20 \\
    --length_required 50 \\
    --detect_adapter_for_pe \\
    --json {log_dir}/{job_tag}_fastp.json \\
    --html {log_dir}/{job_tag}_fastp.html

# Step 2 — Map trimmed reads with Bowtie2
bowtie2 --{args.sensitive} \\
    --phred33 \\
    --threads 8 \\
    --rg-id {ind_name} \\
    --rg PL:ILLUMINA \\
    --rg SM:{ind_name} \\
    -X 500 \\
    -x {bowtie_index} \\
    -1 {path_out}/{ind_name}_trimmed_R1.fastq.gz \\
    -2 {path_out}/{ind_name}_trimmed_R2.fastq.gz | \\
samtools view -bSu - | \\
samtools sort -o {output_bam}
"""
        with open(script_path, 'w') as sh:
            sh.write(slurm_content)

        # Submit the job
        subprocess.run(["sbatch", "--mem=40G", "--cpus-per-task=8", str(script_path)])
        print(f"Submitted job for {ref_name}")

if __name__ == "__main__":
    main()
