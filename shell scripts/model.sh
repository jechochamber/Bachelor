#!/bin/bash
#SBATCH --job-name=model_training_unsampled
#SBATCH --output=/fast/AG_Ohler/jdemoli/out/log_%x.%j.out
#SBATCH --time=72:00:00
#SBATCH --gres=gpu:a40:1
#SBATCH --cpus-per-gpu=8


echo "==================="
date
echo "job: $SLURM_JOB_NAME"
echo "jobID: $SLURM_JOB_ID"
echo "node: $SLURM_JOB_NODELIST"
echo "==================="

# source anaconda (bashrc not evaluated in non-interactive sessions, thus we need to specifically source it)
source ~/.bashrc


# activate desired conda environment
conda activate modelenv
cd '/fast/AG_Ohler/jdemoli/bachelorgit/src'
python3 model_unsampled_data.py