working_dir=$(pwd -P)
slurm_profile="$working_dir/slurm-profile"

export PATH=$PATH:$slurm_profile

calculate_number_of_jobs(){
    declare -i running_jobs
    declare -i available_slots
    running_jobs=$(squeue -u $USER|wc -l)-1
    available_slots=990-$running_jobs
    echo $available_slots
}

execute_cmd(){
    # first arg: dag / dry-tun, / local / cluster
    # second arg: sample folder name
    # third arg: number of jobs
    if [ $1 == "dry-run" ]
    then
        # dry run 
        snakemake -s integrated-pipeline.snake --configfile ./integrated-pipeline-profile/config.yml \
            --verbose -d $working_dir/$2 \
            -np > snakemake-dry-run.log 2>&1
            mv snakemake-dry-run.log $2
    fi

    if [ $1 == "local" ]
    then
        # local computer command
        snakemake -s integrated-pipeline.snake --configfile ./integrated-pipeline-profile/config.yml \
            -d $working_dir/$2 --jobs 1 --use-conda --conda-prefix "$4" > snakemake.log 2>&1
        mv snakemake.log $2
    fi

    if [ $1 == "cluster" ]
    then
        # cluster command
        snakemake -s integrated-pipeline.snake --configfile ./integrated-pipeline-profile/config.yml \
           --verbose -d $working_dir/$2 \
           --jobs $3 --cluster-config ./integrated-pipeline-profile/cluster-config.json  \
           --profile $slurm_profile \
           --immediate-submit --notemp --use-conda --conda-prefix "$4" > snakemake.log 2>&1
            mv snakemake.log $2
    fi
}

samples=($(find sample* -maxdepth 1 -iname "sample*" -type d))

type="cluster"

# if you dont have the conda envs already
# launch 'sbatch ./prepare_envs.sbatch' after adapting conda_prefix variable to your setup
conda_prefix=$(grep 'conda_prefix=' prepare_envs.sbatch|cut -d"'" -f2)

for sample in ${samples[@]}
do
    jobs=$(calculate_number_of_jobs)
    while true
    do
        if [ "$jobs" -lt "10" ] && [ "$type" == "cluster" ]
        then
            sleep 3600
        else
            break
        fi
        jobs=$(calculate_number_of_jobs)
    done
    sed -i -re 's/.+"name":.+/        "name": "'$(basename $PWD)'.'$sample'.{rule}",/' integrated-pipeline-profile/cluster-config.json
    echo "launching $sample"
    execute_cmd $type $sample $jobs $conda_prefix
    echo "$sample done"
done
