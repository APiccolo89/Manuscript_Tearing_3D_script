#!/bin/sh
#SBATCH -J Sync                    # (job name)
#SBATCH -o %j.out                       # (name of the job output file, %j expands to the job name)
#SBATCH -N 1                            # (Number of requested nodes)
#SBATCH --ntasks-per-node 1             # (Number of requested cores per node)
#SBATCH -t 48:00:00                     # (Requested wall time)

python3 Data_Base_Reader.py > log_test.txt
#
