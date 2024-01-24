#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 1
#SBATCH -q regular
#SBATCH -J QSP_Experiment
#SBATCH --mail-user=rr637@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
#SBATCH -A m4138

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
#applications may perform better with --gpu-bind=none instead of --gpu-bind=single:1 
conda activate SULI
# srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=single:1 /global/u2/r/rr637/QSP_codes_archive/uniprot.py
python /global/u2/r/rr637/QSP_codes_archive/search_cpp_demo0.py