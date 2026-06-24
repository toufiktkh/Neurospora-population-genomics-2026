# Neurospora Population Genomics — Internship Scripts (M1 IMHE 2026)

**Author:** Toufik Yahia Takhi  
**Supervisor:** Pierre Gladieux (PHIM, INRAE Montpellier)  
**Institution:** Université de Montpellier — Spécialité : Master 1 Interactions Microorganismes-Hôtes-Environnements (IMHE)  
**Cluster:** IFB (Institut Français de Bioinformatique) SLURM cluster  
**Report title:** *Structure génétique des populations du champignon endophyte pyrophile Neurospora*

---

## Overview

This repository contains all bioinformatic scripts used in the M1 internship report on the population genetic structure of pyrophilous, endophytic *Neurospora* populations sampled from three post-fire sites:

- **Domefire** — Cima Dome, Mojave National Preserve, California, USA (burned Joshua trees, *Yucca brevifolia* var. *jaegeriana*, 2020 Dome Fire)
- **Villeveyrac** — Hérault, France (burned *Cytisus* sp. shrubs, 2021 wildfire)
- **La Grande Motte** — Hérault, France (burned Lauraceae trees, 2021 wildfire)

The analyses characterise species identity, intraspecific population structure, clonal diversity, and nucleotide diversity and divergence across these three sites.

---

## Repository structure

```
internship_scripts/
│
├── 1.phylogenetic_inference_busco/
│   Scripts for BUSCO-based phylogenomic inference using
│   single-copy core orthologues, MAFFT alignment, trimAl
│   trimming, and RAxML-NG maximum-likelihood tree inference.
│
├── 2.Read_mapp_snp_alling_qual_filter/
│   Scripts for read quality trimming (fastp), reference-guided
│   mapping (Bowtie2), SNP/invariant-site calling (bcftools),
│   VCF quality filtering, and QC statistics extraction and
│   visualisation (sequencing depth, mapping quality, QUAL).
│
├── 3.Pop_structure_admixture_analysis/
│   Scripts for PCA (PLINK + ggplot2), Neighbor-Net network
│   inference (SplitsTree4 via vcf2phylip), and sNMF ancestry
│   proportion estimation (LEA R package).
│
├── 4.Clone_id_spatial_analysis/
│   Scripts for clone identification (poppr R package),
│   distance-threshold selection, clone correction, population
│   assignment from sNMF Q-matrices for pixy, and spatial
│   distribution of clonal genotypes across host trees.
│
├── 5.Population_genomic_stat/
│   Scripts for computing and visualising population genomic
│   summary statistics (nucleotide diversity π, absolute
│   divergence Dxy, Tajima's D) using pixy.
│
Generated_figures/
├── 1.Phylogenetic_tree_reference_assemblies/   ML phylogeny (RAxML-NG)
├── 2.Read_mapp_snp_alling_qual_filter/         alignment heatmaps,
│                                                mapping-rate tables (CSV),
│                                                and QC stats plots per site
├── 3.Pop_structure_admixture_analysis/
│   ├── PCA/                                     per-dataset PCA plots
│   ├── SplitsTree/                              Neighbor-Net networks
│   │                                            (.png) + input files
│   │                                            (.nexus, .phy)
│   └── sNMF/                                    best-K barplots, all-K
│                                                barplots, cross-entropy plots
├── 4.Clone_id_spatial_analysis/                clone-correction histograms
│                                                (per site + combined)
└── 5.Population_genomic_stat/                   pixy π, Dxy and Tajima's D plots
```

---

## Software versions

All analyses were run on the IFB SLURM cluster. Modules were loaded via `module load`.

| Software | Version | Use |
|---|---|---|
| fastp | 1.0.1 | Read quality trimming |
| Bowtie2 | 2.3.4.3 | Reference-guided read mapping |
| samtools | 1.9 | BAM sorting and indexing |
| bcftools | 1.16 | SNP/invariant-site calling and VCF filtering |
| htslib / tabix / bgzip | 1.9 | VCF compression and indexing |
| BUSCO | 6.0.0 | Single-copy core gene identification (sordariales_odb12) |
| MAFFT | 7.525 | Multiple sequence alignment |
| trimAl | 1.4.1 | Alignment trimming |
| RAxML-NG | 1.2.2 | Maximum-likelihood phylogenetic inference |
| PLINK | 1.90b6.18 | SNP filtering, LD pruning, PCA |
| vcftools | 0.1.16 | Per-individual missingness statistics |
| vcf2phylip | 2.8 | VCF to PHYLIP/NEXUS conversion for SplitsTree |
| SplitsTree4 | 4.19.2 | Neighbor-Net network visualisation |
| pixy | 2.0.0.beta14 | π, Dxy, Tajima's D estimation |
| R | 4.5.2 | All downstream analysis and visualisation |
| R — vcfR | 1.16.0 | VCF import in R |
| R — poppr | 2.9.8 | Clone identification and correction |
| R — LEA | 1.4.0 | sNMF ancestry estimation |
| R — ggplot2 | 4.0.3 | Visualisation |
| Python | 3 | Alignment processing and batch SLURM script generation |

---

## Pipeline overview

The analyses follow five sequential steps. Each numbered directory corresponds to one step.

```
1. Phylogenetic inference
   BUSCO (genome mode, sordariales_odb12, MetaEuk)
   → per-gene FASTA files
   → MAFFT alignment
   → trimAl trimming (-automated1)
   → filter (gap ≤5%, length ≥2000bp, ≥5 variable sites)
   → partitioned supermatrix
   → RAxML-NG (--all, --bs-trees autoMRE)
   → rooted on *Boothiella tetraspora* in iTOL v7

2. Read mapping, SNP calling and quality filtering
   fastp (Q20, length ≥50bp, --detect_adapter_for_pe)
   → Bowtie2 (--very-sensitive, -X 500)
   → samtools sort/index
   → bcftools mpileup (-d 1000, -Q 20, -a AD,DP) | call (-m, --ploidy 1/2)
   → bcftools filter (FORMAT/DP<3 masked; QUAL<30 or MQ<30 removed)
   → bcftools norm -d all (deduplication)
   → bcftools view --exclude-types indels,mnps,other (SNP extraction)
   → QC stats extraction (bcftools query) and visualisation (R)

3. Population structure analysis
   PLINK (--mind, --geno 0.4, --maf 0.05, --indep-pairwise 50 5 0.5, --pca 10)
   → PCA visualisation (ggplot2)
   → vcf2phylip → SplitsTree4 (Neighbor-Net, uncorrected distances)
   → sNMF (LEA, K=1–10, 10 runs per K, cross-entropy criterion)

4. Clone identification and spatial analysis
   poppr (diss.dist/nLoc, mlg.filter, clonecorrect, strata = ~tree)
   → clone correction (thresholds: DF1=0.03, DF3=0.06, Vill=0.01, LGM=0.02)
   → strict clone distribution analysis (threshold=0.005)
   → population assignment from sNMF Q-matrices (≥90% ancestry threshold)

5. Population genomic summary statistics
   pixy v2.0.0 (--stats pi tajima_d dxy, --window_size 50000)
   → window filter: no_sites ≥ 10000
   → visualisation (ggplot2 violin and sliding-window plots)
```

---

## Datasets

Raw sequencing data (Illumina paired-end, GBS) are not deposited here as they are unpublished. Reference assemblies (E-DF1, E-DF3, Vill-A1-3, neutre-22) are unpublished and were provided by the Gladieux laboratory (PHIM, INRAE Montpellier).

| Dataset | Reference assembly | Species assignment | n isolates (initial) |
|---|---|---|---|
| Domefire E-DF1 | E-DF1 (unpublished) | *Neurospora* sp. | 51 |
| Domefire E-DF3 | E-DF3 (unpublished) | *N. discreta* PS4B | 9 |
| Villeveyrac | Vill-A1-3 | *N. crassa* clade | 43|
| La Grande Motte | neutre-22 | *Neurospora* aff. *tetrasperma* | 20 |

---

## How to reproduce the analyses

1. Clone this repository:
```bash
git clone https://github.com/toufiktkh/neurospora-population-genomics-2026.git
cd neurospora-population-genomics-2026
```

2. Load the required modules on your SLURM cluster (module names may differ):
```bash
module load bowtie2/2.3.4.3 samtools/1.9 fastp/1.0.1 bcftools busco/6.0.0 \
            mafft/7.525 trimal/1.4.1 plink/1.90b6.18 pixy/2.0.0.beta14 r/4.5.2
```

3. Run scripts in the order of the numbered directories (1 → 5). Each directory contains a mix of bash/SLURM scripts (`.sh`) and R scripts (`.R`). SLURM scripts are submitted with `sbatch`; R scripts are run with `Rscript`.

4. Adjust file paths at the top of each script to match your own directory structure before running.

---

## Notes

- Ploidy is set to 1 (haploid) for all datasets except La Grande Motte, which was called with `--ploidy 2` (suspected diploid *N.* aff. *tetrasperma*).
- The all-sites VCF (including invariant positions, produced by `bcftools call -m`) is required as input for pixy. Variants-only VCFs (produced by PLINK filtertaion `-maff 0.05` or extracted with `--exclude-types`) are used for all other analyses.
- Clone correction was applied separately before sNMF and pixy analyses. The clone-corrected isolate lists are written to `*_snps_only_samples_to_keep.txt` files by the clone-correction script.
- The `high_qual_figures/` directory contains the figures used in the internship report at print resolution.

---

## Citation

If you use or adapt these scripts, please cite:

> Takhi T. Y. (2026). *Structure génétique des populations du champignon endophyte pyrophile Neurospora* — M1 IMHE internship report, Université de Montpellier. Supervised by P. Gladieux, PHIM, INRAE Montpellier. Scripts available at: https://github.com/toufiktkh/neurospora-population-genomics-2026

---

## Contact

Toufik-Yahia Takhi — toufikjunior@gmail.com  
Pierre Gladieux (supervisor) — pierre.gladieux@inrae.fr
