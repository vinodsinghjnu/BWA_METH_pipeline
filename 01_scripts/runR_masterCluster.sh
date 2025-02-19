#!/bin/bash
#
#SBATCH --account mutationalscanning
#SBATCH --ntasks-per-node=1
#SBATCH -J NoJobName
#SBATCH --output=xx-%x.%j.out
#SBATCH --error=xx-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vinod.acear@gmail.com

source /home/vinodsingh/miniforge3/bin/activate R.4.2_cloned_mapdemo_env


which R
# Rscript --vanilla $1 -G $2 -C $3 -I $4 -O $5
# R CMD BATCH --vanilla '--args -G $2 -C $3 -I $4 -O $5' $1

#Rscript --vanilla $1 $2 # $1 is rscript and $2 is string of arguments
Rscript --vanilla $@

#sstat  -j $SLURM_JOB_ID.batch --format=JobID,MaxVMSize

#seff $SLURM_JOBID

#sacct --duplicates -j $SLURM_JOBID --format=JobID,JobName,MaxRSS,MaxRSSNode,MaxRSSTask,MaxVMSize,MaxVMSizeNode,MaxVMSizeTask

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch


