#!/bin/bash

# Directories
CWD="/mnt/nfs/home/joao/softwares/pyKVFinder/etc/performance"
PARKVFINDER_INSTALLATION="/mnt/nfs/home/joao/softwares/parKVFinder"

# Nthreads
declare -a nthreads=("1" "2" "4" "8" "12" "16" "20" "24")

################## Benchmarking dataset ##################
printf "[===> Downloading kv1000 dataset\n"
git clone https://github.com/jvsguerra/kv1000.git
mv kv1000/kv1000/pdbs ${CWD}/kv1000
rm -rf kv1000/

################## PARKVFINDER ANALYSIS ##################
printf "[===> parKVFinder benchmarking\n"
mkdir ${CWD}/raw/
mkdir ${CWD}/raw/parKVFinder

for i in "${nthreads[@]}"
do

	printf ">>> ${i} threads\n"

	# Run KVFinder for PDBs
	python ${CWD}/scripts/benchmarking/run_parKVFinder.py ${CWD}/kv1000 ${CWD}/raw/parKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm -r ${CWD}/kv1000/KV_Files/
	rm KVFinder.log

done

################## PYKVFINDER V0.1 ANALYSIS ##################
printf "[===> pyKVFinder v0.1 benchmarking\n"
mkdir ${CWD}/raw/
mkdir ${CWD}/raw/pyKVFinder/

for i in "${nthreads[@]}"
do
	printf ">>> ${i} cores\n"

	# Run KVFinder for PDBs
	python ${CWD}/scripts/benchmarking/run_pyKVFinder.py ${CWD}/kv1000 ${CWD}/raw/pyKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm *.KVFinder.results.toml *.KVFinder.output.pdb *.parameters.toml *.log
done

################## COMPILING BENCHMARKING ##################
printf "[===> Compiling benchmarking\n"
mkdir ${CWD}/data
python ${CWD}/scripts/benchmarking/compile_benchmarking.py ${CWD}/raw ${CWD}/kv1000 ${CWD}/data
