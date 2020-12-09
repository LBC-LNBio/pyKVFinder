#!/bin/bash

# Directories
OUTPUT_DIRECTORY=${1%/}

# Nthreads
declare -a nthreads=("1" "2" "4" "8" "12" "16" "20" "24")

# Require positional command-line arguments
if [ -z "$1" ]; then
    printf "> Missing OUTPUT_DIRECTORY!\n"

	# Help menu
	echo 'Usage:'
	echo -e "\t"'./runbenchmarking.sh OUTPUT_DIRECTORY'

	exit
fi

################## Benchmarking dataset ##################
printf "[===> Downloading kv1000 dataset\n"
git clone https://github.com/jvsguerra/kv1000.git
mv kv1000/kv1000/pdbs ${OUTPUT_DIRECTORY}/data/kv1000
rm -rf kv1000/

################## PARKVFINDER ANALYSIS ##################
printf "[===> parKVFinder benchmarking\n"
mkdir ${OUTPUT_DIRECTORY}/data/raw/
mkdir ${OUTPUT_DIRECTORY}/data/raw/parKVFinder

for i in "${nthreads[@]}"
do

	printf ">>> ${i} threads\n"

	# Run KVFinder for PDBs
	python ${OUTPUT_DIRECTORY}/scripts/benchmarking/run_parKVFinder.py ${OUTPUT_DIRECTORY}/data/kv1000 ${OUTPUT_DIRECTORY}/data/raw/parKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm -r ${OUTPUT_DIRECTORY}/data/kv1000/KV_Files/
	rm KVFinder.log

done

################## PYKVFINDER V0.1 ANALYSIS ##################
printf "[===> pyKVFinder v0.1 benchmarking\n"
mkdir ${OUTPUT_DIRECTORY}/data/raw/pyKVFinder/

for i in "${nthreads[@]}"
do
	printf ">>> ${i} cores\n"

	# Run KVFinder for PDBs
	python ${OUTPUT_DIRECTORY}/scripts/benchmarking/run_pyKVFinder.py ${OUTPUT_DIRECTORY}/data/kv1000 ${OUTPUT_DIRECTORY}/data/raw/pyKVFinder ${i}

	# Save PDBs runs inside KVFiles_ncores
	rm *.KVFinder.results.toml *.KVFinder.output.pdb *.parameters.toml *.log
done

################## COMPILING BENCHMARKING ##################
printf "[===> Compiling benchmarking\n"
python ${OUTPUT_DIRECTORY}/scripts/benchmarking/compile_benchmarking.py ${OUTPUT_DIRECTORY}/data/raw ${OUTPUT_DIRECTORY}/data/kv1000 ${OUTPUT_DIRECTORY}/data
