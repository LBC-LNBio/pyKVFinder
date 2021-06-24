#!/bin/bash

# Running on cluster (example)
# ./compile_parKVFinder.sh ~/softwares/pyKVFinder/tests/performance ~/softwares/parKVFinder

# Directories
OUTPUT_DIRECTORY=${1%/}
PARKVFINDER_LOCATION=${2%/}

# Nthreads
old=("omp_get_num_procs () - 1")
declare -a arr=("1" "2" "4" "8" "12" "16")

if [ -z "$1" ] || [ -z "$2" ]; then

	if [ -z "$1" ]; then
	printf "> Missing OUTPUT_DIRECTORY!\n"
	fi

	if [ -z "$2" ]; then
	printf "> Missing PARKVFINDER_LOCATION!\n"
	fi

	# Help menu
	echo 'Usage:'
	echo -e "\t"'./compile_parKVFinder.sh OUTPUT_DIRECTORY PARKVFINDER_LOCATION'
	exit
fi

# Download parKVFinder
git clone https://github.com/LBC-LNBio/parKVFinder.git
mv parKVFinder ${PARKVFINDER_LOCATION}

for i in "${arr[@]}"
do

	printf ">>> ${i} threads\n"

	# Change KVFinder scripts
	sed -i -e "s/ncores = ${old}/ncores = $i/" ${PARKVFINDER_LOCATION}/src/parKVFinder.c
	sed -i -e "s/ncores = ${old}/ncores = $i/" ${PARKVFINDER_LOCATION}/src/matrixprocessing.c

	# Compile new parKVFinder for i cores
	cd ${PARKVFINDER_LOCATION}; 
	make clean; 
	make; 
	cd ${OUTPUT_DIRECTORY};

    # Prepare executables
    mv ${PARKVFINDER_LOCATION}/parKVFinder ${OUTPUT_DIRECTORY}/scripts/parKVFinder${i}

	# Save ncores as old
	old=($i)

done

# Revert all changes in KVFinder scripts
sed -i -e "s/ncores = ${old}/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_LOCATION}/src/parKVFinder.c
sed -i -e "s/ncores = ${old}/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_LOCATION}/src/matrixprocessing.c

# Compile new KVFinder for i cores
cd ${PARKVFINDER_LOCATION}; make clean; make; cd ${OUTPUT_DIRECTORY}/scripts
