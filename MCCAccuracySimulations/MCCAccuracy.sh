#!/bin/sh
#Reserve 2 CPUs for this job
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Request it to run this for HH:MM:SS with ?G per core
#SBATCH --time=05:59:00
# Specify the queue (30min, 6hours, 1day, 1week)
#SBATCH --qos=6hours
#
#SBATCH --output=SLURM_%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=SLURM_%j.err                  # where to store error messages

source ~/miniconda3/etc/profile.d/conda.sh
conda activate MTK_flu
export PATH="$HOME/.julia/bin:$PATH"

# execute whatever was passed to the script
snakemake --unlock
snakemake -j 
