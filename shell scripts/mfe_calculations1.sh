#!/bin/bash
#SBATCH --job-name=MFE_160GC50
#SBATCH --output=/fast/AG_Ohler/jdemoli/out/log_%x.%j.out
#SBATCH --time=48:00:00
#SBATCH --mem=25G

echo "==================="
date
echo "job: $SLURM_JOB_NAME"
echo "jobID: $SLURM_JOB_ID"
echo "node: $SLURM_JOB_NODELIST"
echo "==================="

# source anaconda (bashrc not evaluated in non-interactive sessions, thus we need to specifically source it)
source ~/.bashrc


# activate desired conda environment
conda activate testenv
cd /fast/AG_Ohler/jdemoli/bachelorgit/
python3 mfe_calculations_single.py