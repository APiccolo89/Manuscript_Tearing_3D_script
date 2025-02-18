#!/bin/bash
#SBATCH --job-name Create_Markers_LaMEM_run_mode # (job name)
#SBATCH --output job.%j.out # (name of the job output file, %j expands to the job name)
#SBATCH --ntasks 128 # (Number of requested cores) ! at all
#SBATCH --partition=WAP
#SBATCH --exclude=r12g03 # couldn't find panfs

#module load inteloneapi/mpi/2021.5.0 petsc/3.16.4

mpirun --map-by socket -np $SLURM_NTASKS ../lamem_new/bin/opt/LaMEM -ParamFile SSB_Nevena.dat   > This.txt
#mpirun -n 128 ../lamem/bin/opt/LaMEM -ParamFile SSB_Nevena.dat>job.txt

