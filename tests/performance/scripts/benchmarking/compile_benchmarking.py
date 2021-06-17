import os
import sys
import pandas as pd

# sys.argv[1]: raw directory
# sys.argv[2]: kv1000 directory
# sys.argv[3]: data/ output directory

def get_number_of_atoms(pdb):
    from Bio.PDB import PDBParser
    # Read pdb
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(f"{pdb.replace('kv1000/', '').replace('.pdb', '')}", pdb)
    # Count number of atoms
    n_atoms = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    n_atoms += 1
    return n_atoms

def process_raw_data(fn = sys.argv[1], kv1000 = sys.argv[2], output_directory = sys.argv[3]):   
    # Empty dataframe
    data = {}
    # Get pdb names
    data['pdb'] = [pdb.replace('.pdb', '') for pdb in sorted(os.listdir(kv1000)) if pdb.endswith('.pdb')]
    # Get number of atoms
    data['natoms'] = [get_number_of_atoms(f"{kv1000}/{pdb}.pdb") for pdb in data['pdb']]
    # Get results directories
    d = [d for d in os.listdir(fn) if os.path.isdir(os.path.join(fn, d))]
    # Get time for each software (pyKVFinder and parKVFinder) for different number of threads
    for software in sorted(d):
        for nthreads in sorted(os.listdir(os.path.join(fn, software))):
            data[f"{software}_{nthreads.replace('.csv', '')}"] = pd.read_csv(os.path.join(fn, software, nthreads))[' avg_time']
    # To pandas DataFrame
    data = pd.DataFrame(data)
    data.to_csv(f"{output_directory}/time.csv")
    
    return len(d)

if __name__ == "__main__":
    nsoftware = process_raw_data()
    os.system(f"Rscript {os.path.abspath(os.path.dirname(__file__))}/plot_benchmarking.R {sys.argv[3]} {nsoftware}")