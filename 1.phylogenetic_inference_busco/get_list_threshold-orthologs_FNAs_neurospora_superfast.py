import os
import argparse
import shutil
from multiprocessing import Pool
from collections import defaultdict
from tqdm import tqdm  # External library for progress bars

def process_isolate(isolate_data):
    # This runs inside the worker processes
    isolate, base_path, out_path = isolate_data
    isolate_dir = os.path.join(out_path, isolate)

    os.makedirs(isolate_dir, exist_ok=True)

    # Path construction
    mypath = os.path.join(base_path, isolate)

    found_genes = []
    if os.path.exists(mypath):
        for file in os.listdir(mypath):
            if file.endswith('.fna'):
                gene = file.split('.')[0]
                found_genes.append(gene)
                shutil.copy2(os.path.join(mypath, file), isolate_dir)

    return isolate, found_genes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True)
    parser.add_argument("--included", required=True)
    parser.add_argument("--pathout", required=True)
    parser.add_argument("--cutoff", type=float, default=1.0)
    parser.add_argument("--SCOlist", required=True)
    args = parser.parse_args()

    # 1. Read isolates
    with open(args.included, 'r') as f:
        isolates = [line.strip() for line in f if line.strip()]

    # 2. Setup output directory
    if os.path.exists(args.pathout):
        print(f"Cleaning output directory: {args.pathout}")
        shutil.rmtree(args.pathout)
    os.makedirs(args.pathout)

    # 3. Parallel Processing with Progress Bar
    worker_data = [(iso, args.path, args.pathout) for iso in isolates]

    sco_map = defaultdict(list)
    results = []

    print(f"Processing {len(isolates)} isolates using parallel workers...")

    # Using 'imap_unordered' with tqdm for real-time progress updates
    with Pool() as pool:
        # total=len(isolates) allows tqdm to calculate the % and ETA
        for res in tqdm(pool.imap_unordered(process_isolate, worker_data), total=len(isolates), desc="Copying Genes"):
            results.append(res)

    # 4. Aggregate Results
    print("Filtering Single Copy Orthologs (SCOs)...")
    for isolate, genes in results:
        for gene in genes:
            sco_map[gene].append(isolate)

    # 5. Filter and Write SCOlist
    cutoff_val = args.cutoff * len(isolates)
    with open(args.SCOlist, 'w') as out:
        # Find genes that meet the cutoff
        valid_genes = [gene for gene, isos in sco_map.items() if len(isos) >= cutoff_val]

        # tqdm again if the list is huge, otherwise simple print is fine
        for gene in tqdm(valid_genes, desc="Writing SCO List"):
            out.write(f"{gene}\n")

    print(f"Done! Found {len(valid_genes)} genes meeting the {args.cutoff} cutoff.")

if __name__ == "__main__":
    main()
