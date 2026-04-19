#!/bin/bash
#SBATCH --job-name=rnaseq_fuzzy_cnar
#SBATCH --mail-user=[your email]
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --partition=total
#SBATCH --mem 20G
##SBATCH --time=0-00:00:00
#SBATCH --cpus-per-task=125
##SBATCH --nodelist=stat01
##SBATCH --hint=nomultithread
##SBATCH --output=user_job%j.out
##SBATCH --error=user_error_%j.err

CONTAINER= ../mySing.sif

##export OMP_NUM_THREADS=1
##export OPENBLAS_NUM_THREADS=1
##export MKL_NUM_THREADS=1

srun -c $SLURM_CPUS_PER_TASK /usr/bin/singularity exec ${CONTAINER} R CMD BATCH rnaseq_fuzzy/run_sim.R
