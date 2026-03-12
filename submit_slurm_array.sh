#!/bin/bash
#SBATCH --job-name=md_grid_sweep
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-23          # 24 jobs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Slurm

# Parameter grid (6 temps x 4 angles = 24 combinations)
TEMPS=(0.6 0.8 1.0 1.2 1.4 1.6)
ANGLES=(5.0 10.0 15.0 20.0)

# Map Slurm Array ID to correct parameter combination
# Division gets the Temp index (0,0,0,0,0,0, 1,1,1... etc)
TEMP_IDX=$(( SLURM_ARRAY_TASK_ID / 4 ))
# Modulo gets the Angle index (0,1,2,3,4,5, 0,1,2... etc)
ANGLE_IDX=$(( SLURM_ARRAY_TASK_ID % 4 ))

CURRENT_TEMP=${TEMPS[$TEMP_IDX]}
CURRENT_ANGLE=${ANGLES[$ANGLE_IDX]}

# Workspace setup
BASE_DIR=$SLURM_SUBMIT_DIR
JOB_WORKSPACE="${BASE_DIR}/runs/job_${SLURM_ARRAY_JOB_ID}/task_${SLURM_ARRAY_TASK_ID}"

mkdir -p ${JOB_WORKSPACE}

rsync -a --exclude 'build' \
         --exclude 'results' \
         --exclude '.venv' \
         --exclude 'runs' \
         --exclude '*.sif' \
         ${BASE_DIR}/ ${JOB_WORKSPACE}/

cd ${JOB_WORKSPACE}
#ln -s ${BASE_DIR}/.venv ${JOB_WORKSPACE}/.venv

# Modify input.dat (Line 2 for Temp, Line 11 for maxDih)
# We use '|' as the sed delimiter so the slashes in the comments don't break it
sed -i "2s|.*|${CURRENT_TEMP}                     ! temperature (in units of epsilon/k_B; k_B = 0.0019872041 kcal/mol*K)|" input.dat
sed -i "11s|.*|${CURRENT_ANGLE}                     ! maxDih: maximum dihedral angle variation (degrees)|" input.dat

# Run simulation
apptainer exec ${BASE_DIR}/md_tools.sif make MODE=sequential plot