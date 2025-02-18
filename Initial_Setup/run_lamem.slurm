#!/bin/bash
#SBATCH -J Create_Markers_LaMEM_run_mode # (job name)
#SBATCH -o job.%j.out # (name of the job output file, %j expands to the job name)
#SBATCH -N 8 # (Number of requested nodes)
#SBATCH -n 128 # (Number of requested cores per node)
#SBATCH -t 24:00:00 # (Requested wall time)

module load inteloneapi/mpi/2021.5.0
module load petsc/3.16.4

mpirun -n 128 ../lamem/bin/opt/LaMEM -ParamFile SSB_Nevena.dat  >job_restart.txt 
