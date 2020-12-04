#!/bin/bash

# Directories
CWD="/mnt/nfs/home/joao/softwares/pyKVFinder/etc/performance"
PARKVFINDER_INSTALLATION="/mnt/nfs/home/joao/softwares/parKVFinder"

# Nthreads
old=("omp_get_num_procs () - 1")
declare -a arr=("1" "2" "4" "8" "12" "16" "20" "24")

for i in "${arr[@]}"
do

	printf ">>> ${i} threads\n"

	# Change KVFinder scripts
	sed -i -e "s/ncores = ${old}/ncores = $i/" ${PARKVFINDER_INSTALLATION}/src/parKVFinder.c
	sed -i -e "s/ncores = ${old}/ncores = $i/" ${PARKVFINDER_INSTALLATION}/src/matrixprocessing.c

	# Compile new parKVFinder for i cores
	cd ${PARKVFINDER_INSTALLATION}; 
	make clean; 
	make; 
	cd ${CWD};

    # Prepare executables
    mv ${PARKVFINDER_INSTALLATION}/parKVFinder ${PARKVFINDER_INSTALLATION}/parKVFinder${i}

	# Save ncores as old
	old=($i)
done

# Revert all changes in KVFinder scripts
sed -i -e "s/ncores = 24/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_INSTALLATION}/src/parKVFinder.c
sed -i -e "s/ncores = 24/ncores = omp_get_num_procs () - 1/" ${PARKVFINDER_INSTALLATION}/src/matrixprocessing.c

# Compile new KVFinder for i cores
cd ${PARKVFINDER_INSTALLATION}; make clean; make; cd ${CWD}