**COD-dipp: Closed Open De novo – deep immunopeptidomics pipeline, An immunopeptidomics integrated pipeline**

# Introduction
This pipeline has been developed as an effort to characterize the MHC system across 26 already published immunopeptidomic datasets. In more details,the approach relies on 3 robust mass spectrometry identification strategies (closed search, open search and denovo). Thus, allowing the detection of the wild type peptides (closed search), the post translational landscape (open search) and novel sequences (denovo) even when complementary transcriptomic/genomic data is not available.

# How to use

## Step 1: prerequisite
 - Ensure you have anaconda installed <https://www.anaconda.com/distribution>
 - Ensure you have snakemake installed (python3 -m pip install --user 'snakemake==5.4.5')
 - Download msfragger from this link: <http://msfragger.nesvilab.org> and update its absolute path in integrated-pipeline-profile/config.yml (V2.2 has been tested only)
 - Create a Target-decoy database with a decoy tag "rev_" (only ENSEMBL DATABASE HEADER ">ENSP...|ENST..|ENSG... DESCRIPTION" has been tested)
 - Generate the resource files needed following the [resources generation tutorial](./scripts/prepare_annotation) or download from figshare
 - Edit the [integrated-pipeline-profile/config.yml](./integrated-pipeline-profile/config.yml) configuration file
 - Edit [integrated-pipeline-profile/cluster-config.json](./integrated-pipeline-profile/cluster-config.json) with your own Account name on a slurm cluster and the resources compatible with the HPC settings

## Step 2: Preparing the environment
```
# create a directory for your study
study_id="Example_Study"
mkdir $study_id && cd $study_id

# clone environement
git clone ssh://git@git.plgrid.pl:7999/iccvs/ipip.git
cd ipip
rm -rf sample_test database.fasta
mv $(ls -A) ../
cd ..
rmdir ipip

# copy files
unzip pipeline_annotation_files.zip # link: <https://doi.org/10.6084/m9.figshare.16538097>
resource_files_dir="/PATH/TO/pipeline_annotation_files/folder"
# copy your database to your study working dir
cp $resource_files_dir/2019-04-30-td-Homo_sapiens_GRCh38_biomart.fasta ./database.fasta
# copy knapsack.npy (optional to speed up denovo, if not available DeepNovoV2 will generate it)
cp -lr $resource_files_dir/knapsack.npy scripts/DeepNovoV2
```

## Step 3: install conda environements
 - edit "conda_prefix" path in prepare_envs.sbatch to your desired location on the cluster
 - edit the ./prepare_envs.sbatch SLURM "-A" parameter
 - launch this command to create the conda environement at the "conda_prefix" location using the command below
```
sbatch ./prepare_envs.sbatch
```

## Step 4: prepare the samples
Make sure the samples are organized as the example below

```
working_dir
    │
    └── sample_nameofsample1
        ├── denovo
        │   ├── fraction1.mgf
        │   └── fraction2.mgf
        ├── fraction1.mzML
        └── fraction2.mzML
```

NOTE:
- Ensure using MSConvert proteowizard version >= 3.0.19304.503cb4044 for mgf and mzML.
- for MSconvertGUI make sure the "TPP compatibility" is checked and PeakPicking is enabled
- for the MSconvert command line version, please take a look at the commands below

```
msconvert fraction1.raw --mgf --filter "peakPicking true 1-" --filter 'titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:"<SourcePath>", NativeID:"<Id>"'
msconvert fraction1.raw --mzML --filter "peakPicking true 1-" --filter 'titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:"<SourcePath>", NativeID:"<Id>"'
```

## Step 5: launch the pipeline
There are 3 modes to launch the pipeline:
 - dry-run mode: will only print the rules and commands that will be run in a log file and will store it inside the sample folder under snakemake-dry-run.log
 - local mode: will run the pipeline on the local machine. This is useful for testing or running a specific rule inside an interactive job
 - cluster mode: will use the SLURM workload manager to launch the jobs in the background
```
# the following scripts assume that the samples folder names start with 'sample_'
# you could edit the type variable in launch_pipeline.sh to either "cluster", "local", "dry-run"
bash launch_pipeline.sh
```

# reports and logs
Execution reports
 - The pipeline will generate the following logs:
   - slurm-logs director containing all the launched slurm jobs (cmd, out and err files per job)
   - snakemake-dry-run.log *verbose log without executing jobs*
   - snakemake.log *the job execution log*

```
working_dir
    │
    └── sample/slurm-logs
        ├── rulename-date-time_cmd.log
        ├── rulename-date-time_err.log
        ├── rulename-date-time_out.log
        ├── snakemake-dry-run.log
        └── snakemake.log
```

Please check scripts/reproting for concatenation, reporting and plotting script

# Authors:
 - Georges Bedran: gbadran_90@live.com
