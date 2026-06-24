import os
import argparse
from multiprocessing import Pool
from collections import defaultdict
from tqdm import tqdm

def read_gene_file(file_info):
    """Worker function to read a single fasta file and return the formatted content."""
    dir_name, file_name, path_to_dir, scos_set = file_info
    gene = file_name.split('.')[0]
    
    if gene not in scos_set:
        return None
    
    content = []
    file_path = os.path.join(path_to_dir, file_name)
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                content.append(f">{dir_name}\n")
            else:
                content.append(line)
    return gene, "".join(content)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True)
    parser.add_argument("--pathout", required=True)
    parser.add_argument("--genelist", required=True)
    args = parser.parse_args()

    os.makedirs(args.pathout, exist_ok=True)

    # 1. Use a SET for O(1) lookup speed (massive improvement over lists)
    with open(args.genelist, 'r') as f:
        scos_set = {line.strip() for line in f if line.strip()}

    # 2. Prepare task list for parallel workers
    print("Scanning directories...")
    tasks = []
    subdirs = [d for d in os.listdir(args.path) if os.path.isdir(os.path.join(args.path, d))]
    
    for dir_ in subdirs:
        dir_path = os.path.join(args.path, dir_)
        for file in os.listdir(dir_path):
            tasks.append((dir_, file, dir_path, scos_set))

    # 3. Parallel Read
    # This saturates your CPU cores to bypass the I/O bottleneck
    sequences = defaultdict(list)
    print(f"Processing {len(tasks)} files across {os.cpu_count()} cores...")
    
    with Pool() as pool:
        # imap_unordered is the fastest way to stream results back
        for result in tqdm(pool.imap_unordered(read_gene_file, tasks), total=len(tasks), desc="Reading"):
            if result:
                gene, data = result
                sequences[gene].append(data)

    # 4. Fast Write
    print(f"Writing {len(sequences)} gene fasta files...")
    for gene, contents in tqdm(sequences.items(), desc="Writing"):
        out_path = os.path.join(args.pathout, f"{gene}.fasta")
        with open(out_path, 'w') as f:
            f.write("".join(contents))

if __name__ == "__main__":
    main()
