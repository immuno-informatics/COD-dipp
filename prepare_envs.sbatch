#!/bin/bash
#SBATCH -A AccountID
#SBATCH -n 1
#SBATCH --mem 10GB
#SBATCH --time=01:00:00
#SBATCH --output=prepare_envs.%J.out

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1


# And finally run the jobs
conda_prefix='PATH/TO/DESIRED/CONDA/ENV/OTUPUT/DIR'

cp -r ./scripts/Scavager-0.1.9_combined $conda_prefix
echo "preparing envs"
snakemake -s integrated-pipeline.snake --configfile ./integrated-pipeline-profile/config.yml --create-envs-only --use-conda --conda-prefix $conda_prefix
echo "preparing envs done"
rm -rf .snakemake
rm -rf $conda_prefix/Scavager-0.1.9_combined
echo "env created at $conda_prefix"
