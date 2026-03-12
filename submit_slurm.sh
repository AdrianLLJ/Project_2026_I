#!/bin/bash
#SBATCH --job-name=md_local_run
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# Get CWD
BASE_DIR=$SLURM_SUBMIT_DIR
JOB_WORKSPACE="${BASE_DIR}/runs/job_${SLURM_JOB_ID}"

# Create run workspace
mkdir -p ${JOB_WORKSPACE}

rsync -a --exclude 'build' \
         --exclude 'results' \
         --exclude '.venv' \
         --exclude 'runs' \
		 --exclude '*.sif' \
         ${BASE_DIR}/ ${JOB_WORKSPACE}/

cd ${JOB_WORKSPACE}
ln -s ${BASE_DIR}/.venv ${JOB_WORKSPACE}/.venv

apptainer exec ${BASE_DIR}/md_tools.sif make clean
apptainer exec ${BASE_DIR}/md_tools.sif make MODE=sequential plot