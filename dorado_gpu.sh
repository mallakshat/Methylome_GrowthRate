#!/bin/bash
#SBATCH --job-name=dorado_gpu_job        # Job name
#SBATCH --error=dorado_gpu_error.txt     # error file
#SBATCH --ntasks=1                       # Number of tasks (single task)
#SBATCH --nodes=1       	             # Request 1 nodes
#SBATCH --gres=gpu:2    	             # 2 GPU per node
#SBATCH --partition=gpu-8                # Use the 'gpu-8' partition
#SBATCH --time=24:00:00                  # Time limit (1 day)

module load cuda
export CUDA_VISIBLE_DEVICES=$(echo $SLURM_JOB_GPUS | tr ',' ' ')

# Run Dorado with GPU
~/dorado-1.3.1-linux-x64/bin/dorado basecaller ~/dorado-1.3.1-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
	--device cuda:all \
    sample.pod5 \
    --modified-bases-models ~/dorado-1.3.1-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1  \
    > sample.bam
