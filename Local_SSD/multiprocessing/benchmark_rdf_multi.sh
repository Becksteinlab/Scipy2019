#!/bin/bash
#SBATCH -J PMDA_BM                  # name
#SBATCH --partition=compute
#SBATCH --nodes=1                        # Total number of nodes requested (16 cores/node). You may delete this line if wanted
#SBATCH --ntasks-per-node=24            # Total number of mpi tasks requested
#SBATCH --export=ALL
#SBATCH -t 08:00:00                      # wall time (D-HH:MM)
#SBATCH --mail-type=ALL                # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=sfan19@asu.edu  # send-to address

bash /home/sfan19/.bashrc

echo $SLURM_JOB_ID
echo $USER

python benchmark_rdf_multi.py /scratch/$USER/$SLURM_JOB_ID
