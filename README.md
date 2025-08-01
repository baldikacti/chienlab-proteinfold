# Introduction

This is a convienience [Nextflow](https://www.nextflow.io/) pipeline for generating high-throughtput bait:pair protein structure predictions using [Colabfold](https://github.com/sokrypton/ColabFold), or [Alphafold3](https://github.com/google-deepmind/alphafold3). The pipeline simply accepts a `tsv` file of bait and prey sequences as described in the `Usage examples` section and produces structure predictions and a ranked result file using the `ipTM` values. 

This pipeline is in active development. All versions will be tagged with a version number, so you can use the `nextflow run -r <version>` to use specific versions of the pipeline.

This pipeline is meant for the [Unity HPC](https://unity.rc.umass.edu/index.php) cluster users using [Slurm](https://slurm.schedmd.com/documentation.html) scheduler. It can be adjusted to run locally or other schedulers relatively easily. The pipeline is inspired by several other open source projects like [nf-core/proteinfold](https://github.com/nf-core/proteinfold), and[AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown).

# Usage examples

## Usage: ColabFold

1. Download Proteome `tsv` file for your organism of interest from Uniprot. (Optional)
    - Navigate to `https://www.uniprot.org/`
    - Select `Proteomes` from the the dropdown menu from left side of the search bar
    - Enter the name of your organism of interest
    - From the results, click on the hyperlink under the **Entry** tab. (Ex. **UP000005640** for homo spiens)
    - Scroll down to **Components** section
    - Click `Download` top left of this section
    - Select `Download only reviewed` or `Download reviewed and unreviewed`
    - Select **Format** = `TSV`
    - Select **Compressed** = `No`
    - Optional: You can customize columns in the data by adding or substracting whatever you want. This will be automatically added the the `ranked_results` table at the end of the pipeline

2. Create a new `tsv` file with `Entry` and `bait` columns as below. `Entry` column can be either a UniprotID or a relative path to a `fasta` file. `fasta` file/s can have multiple entries.

Set `bait` = 1 for your bait protein/s. And 0 for every pair you want generated. 

**acclist.tsv**
| Entry     | bait  |
| :-------: | :--:  |
| P41797    | 1     |
| P25685    | 0     |
| Q9UNE7    | 0     |
| Q02790    | 0     |

**OR**

**acclist.tsv**
| Entry            | bait  |
| :--------------: | :--:  |
| ref/P41797.fasta | 1     |
| P25685           | 0     |
| Q9UNE7           | 0     |
| Q02790           | 0     |

3. Example submission script below.

**Mandatory arguments:**

- **input** = Path to `tsv` file. [`/path/to/input.tsv`]

- **output** = Path to result directory. [`/path/to/results`]

- **mode** = Sets the prediction mode. [`colabfold`]

- **num_recycles_colabfold** = Number of recycles to use in Colabfold. Higher the number better the prediction, but the slower the pipeline. [integer]


**Optional arguments:**

- **org_ref** = Path to the `tsv` file downloaded from Uniprot. [`/path/to/uniprot_organism_reference.tsv`]

- **top_rank** = Number of top ranked (by `ipTM`) `bait:pair` predictions to pick for rerunning with 20 recycles for better prediction quality. [integer]

- **APPTAINER_CACHEDIR** = Absolute path to where `apptainer` will cache the containers used in the pipeline. [`/path/to/.apptainer/cache`]

- Additional optional paramaters can be found in `nextflow.config` file.


**Nextflow arguments:**

- **profile** = Set the dependency management profile. [Institution/docker/singularity/apptainer]

- **resume** = Enables the pipeline to be used repeatedly. The pipeline will only run incomplete processes when rerun with the same inputs.


**main.sh**
```bash
#!/usr/bin/bash
#SBATCH --job-name=chienlab-proteinfold     # Job name
#SBATCH --partition=cpu                     # Partition (queue) name
#SBATCH -c 2                                # Number of CPUs
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --mem=10gb                          # Job memory request
#SBATCH --time=14-00:00:00                  # Time limit days-hrs:min:sec
#SBATCH -q long
#SBATCH --output=logs/chienlab-proteinfold_%j.log

module load nextflow/24.10.3 apptainer/latest

APPTAINER_CACHEDIR=/path/to/.apptainer/cache  # Path to cache directory for apptainer cache

export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

nextflow run baldikacti/chienlab-proteinfold -r v0.7.4 \
      --input /path/to/acclist.tsv \
      --outdir /path/to/results \
      --org_ref /path/to/organism_reference.tsv \
      --mode colabfold \
      --num_recycles_colabfold 5 \ 
      --top_rank 2 \
      --colabfold_model_preset "alphafold2_multimer_v3" \
      -profile unity \
      -resume
```

4. Submit to slurm with `sbatch main.sh`

## Results

### Directory Structure

Example results output directory structure can be found below. 

**Legend**

- *colabfold*: Contains the `ColabFold` prediction results for each `bait:prey` pair.

- *screen*: Contains the prediction pairs based on the initial screen using the recycle count from `num_recycles_colabfold` paramater

- *toprank*: Contains the prediction pairs from top ranked pairs based on `ipTM` score with 20 recycles. `top_rank` flag sets how many pairs should be rerun with 20 recycles.

- *fasta*: Contains all `bait:prey` combined FASTA files.

- *pipeline_info*: Contains pipeline execution summaries

- *ranked_results*: Contains the `ranked_results.tsv` file that contains ranked (by `ipTM`) `bait:prey` predictions and optional protein information.

```bash
<outdir>
├── colabfold
│   ├── screen
│   │   ├── # Directories for each bait:prey pair with ColabFold results
│   │   ...
│   └── toprank
│       ├── # Directories for each bait:prey pair with ColabFold results from top ranked pairs
│       ...
├── fasta
│   ├── # Paired FASTA files
│   ...
├── pipeline_info
│   ├── # Execution summaries (html, txt)
│   ...
└── ranked_results
    └── ranked_results.tsv
```

## Usage: Alphafold3

1. Create a new `tsv` file with `Entry` and `bait` columns as below. `Entry` column can be a UniprotID, a relative path to a `fasta` file, a `CCD` code prefixed by `CCD:`, or a `SMILES` code prefixed by `SMILES:`.

Set `bait` = 1 for your bait protein/s. And 0 for every pair you want generated. 

**acclist.tsv**
| Entry         | bait  |
| :-----------: | :--:  |
| CCD:MG        | 1     |
| P25685        | 0     |
| ref/my.fasta  | 0     |
| Q02790        | 0     |

2. Example submission script below.

**Mandatory arguments:**

- **input** = Path to `tsv` file. [`/path/to/input.tsv`]

- **output** = Path to result directory. [`/path/to/results`]

- **model_dir** = Path to the Alphafol3 model paramters directory. [`/path/to/params`]
    
- **db_dir** = Path to the Alphafold3 database. [`/datasets/bio/alphafold3`]

- **mode** = Sets the prediction mode. [`alphafold3`]


**Optional arguments:**

- **max_template_date** = The date for max template date to be used. [`2021-09-30`]

- **num_recycles** = Number of recycles to be used for inference. [10]

- **inf_batch** = Number used for batching number of inference runs per GPU. Used for efficiency. [20]

- **APPTAINER_CACHEDIR** = Absolute path to where `apptainer` will cache the containers used in the pipeline. [`/path/to/.apptainer/cache`]

- Additional optional paramaters can be found in `nextflow.config` file.


**Nextflow arguments:**

- **profile** = Set the dependency management profile. [Institution/docker/singularity/apptainer]

- **resume** = Enables the pipeline to be used repeatedly. The pipeline will only run incomplete processes when rerun with the same inputs.


**main.sh**
```bash
#!/usr/bin/bash
#SBATCH --job-name=chienlab-proteinfold     # Job name
#SBATCH --partition=cpu                     # Partition (queue) name
#SBATCH -c 2                                # Number of CPUs
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --mem=10gb                          # Job memory request
#SBATCH --time=14-00:00:00                  # Time limit days-hrs:min:sec
#SBATCH -q long
#SBATCH --output=logs/chienlab-proteinfold_%j.log

module load nextflow/24.10.3 apptainer/latest

APPTAINER_CACHEDIR=/path/to/.apptainer/cache  # Path to cache directory for apptainer cache

export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

nextflow run baldikacti/chienlab-proteinfold -r v0.7.4 \
      --input /path/to/acclist.tsv \
      --outdir /path/to/results \
      --model_dir /path/to/model/params \
      --db_dir /path/to/db \
      --mode alphafold3 \
      -profile unity \
      -resume
```

3. Submit to slurm with `sbatch main.sh`

# Results

## Directory Structure

Example results output directory structure can be found below. 

**Legend**

- *alphafold3*: Contains the `Alphafold3` prediction results for each `bait:prey` pair.

- *input_json*: Contains all `bait:prey` combined JSON files structured in Alphafold3 required structure.

- *msa*: Contains the MSAs generated from input JSON files

- *folds*: Contains directories for each inference result

- *pipeline_info*: Contains pipeline execution summaries

- *ranked_results*: Contains the `ranked_results.tsv` file that contains ranked (by `ipTM`) `bait:prey` predictions.

```bash
<outdir>
├── colabfold
│   ├── input_json # Contains the generated input JSON files for alphafold3
│   │   ...
│   ├── msa # Contains the MSAs generated from input JSON files
│   │   ...
│   └── folds # Contains directories for each inference result
│       ...
├── pipeline_info
│   ├── # Execution summaries (html, txt)
│   ...
└── ranked_results
    └── ranked_results.tsv
```

## Pipeline Summary

When a run successfully finishes, the `.log` file (set by `#SBATCH --output=/path/to/mylog_%j.log`) will contain a short summary of total execution time, successful and failed jobs. (Check `pipeline_info` directory for detailed execution summaries.)

The pipeline is set up in a way if a job fails, it is retried with higher resources up to two times and `ignored` on the 3 failed attempt. Due to this, if you see a number next to the `Ignored` section in the summary it means that that number of predictions completely failed to produce results. This is most likely due to the size limitation of `ColabFold` or `Alphafol3`.

```bash
Completed at: <Date and time when the pipeline finished>
Duration    : <Total execution duration>
CPU hours   : <CPU hours>
Succeeded   : <number of successfull jobs>
Failed      : <number of failed jobs>
Ignored     : <number of ignored jobs>
```
