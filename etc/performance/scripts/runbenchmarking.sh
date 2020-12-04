#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=pyKVFinder
#SBATCH --output=serial.log
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=24
#SBATCH -w cpu001

# Directories
CWD="/mnt/nfs/home/joao/softwares/pyKVFinder/etc/performance"
PARKVFINDER_INSTALLATION="/mnt/nfs/home/joao/softwares/parKVFinder"

# Nthreads
old=("omp_get_num_procs () - 1")
declare -a arr=("1" "2" "4" "8" "12" "16" "20" "24")

################## Benchmarking dataset ##################
printf "[===> Downloading kv1000 dataset\n"
git clone https://github.com/jvsguerra/kv1000.git
mv kv1000/kv1000/pdbs ${CWD}/kv1000
rm -rf kv1000/

################## PARKVFINDER ANALYSIS ##################
printf "[===> parKVFinder benchmarking\n"
mkdir ${CWD}/raw/
mkdir ${CWD}/raw/parKVFinder

for i in "${arr[@]}"
do

	printf ">>> ${i} threads\n"

	# Run KVFinder for PDBs
	python ${CWD}/scripts/run_parKVFinder.py ${CWD}/kv1000 ${CWD}/raw/parKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm -r kv1000/KV_Files/

	# Save ncores as old
	old=($i)
done

# Revert all changes in KVFinder scripts
sed -i -e "s/ncores = 24/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_INSTALLATION}/src/parKVFinder.c
sed -i -e "s/ncores = 24/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_INSTALLATION}/src/matrixprocessing.c

# Compile new KVFinder for i cores
cd ${PARKVFINDER_INSTALLATION}; make clean; make; cd ${CWD}

################## PYKVFINDER V0.1 ANALYSIS ##################
printf "[===> pyKVFinder v0.1 benchmarking\n"
mkdir ${CWD}/raw/
mkdir ${CWD}/raw/pyKVFinder/

for i in "${arr[@]}"
do
	printf ">>> ${i} cores\n"

	# Run KVFinder for PDBs
	python ${CWD}/scripts/run_pyKVFinder.py ${CWD}/kv1000 ${CWD}/raw/pyKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm *.KVFinder.results.toml *.KVFinder.output.pdb *.parameters.toml *.log
done

################## COMPILING BENCHMARKING ##################
printf "[===> Compiling benchmarking\n"