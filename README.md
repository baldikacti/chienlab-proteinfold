# Introduction

This is a convienience [Nextflow](https://www.nextflow.io/) pipeline for generating high-throughtput bait:pair protein structure predictions using [Colabfold](https://github.com/sokrypton/ColabFold). The pipeline simply accepts a `tsv` file of UniprotIDs and bait status as described in the `Usage` section and produces structure predictions and a ranked result file using the `ipTM` values. 

This pipeline is in active development. All versions will be tagged with a version number, so you can use the `nextflow run -r <version>` to use specific versions of the pipeline.

This pipeline is meant for the [Unity HPC](https://unity.rc.umass.edu/index.php) cluster users using [Slurm](https://slurm.schedmd.com/documentation.html) scheduler. It can be adjusted to run locally or other schedulers relatively easily. The pipeline is inspired by several other open source projects like [nf-core/proteinfold](https://github.com/nf-core/proteinfold), and [LazyAF](https://github.com/ThomasCMcLean/LazyAF).

# Usage

1. Download Proteome `tsv` file for your organism of interest from Uniprot.
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

2. Create a new `tsv` file with `Entry` and `bait` columns as below. `Entry` column needs to be either the full list or a subset of the `Entry` column in the proteome `tsv` file downloaded in step 1. 

Set `bait` = 1 for your bait protein/s. And 0 for every pair you want generated. 

**acclist.tsv**
| Entry     | bait  |
| :-------: | :--:  |
| P41797    | 1     |
| P25685    | 0     |
| Q9UNE7    | 0     |
| Q02790    | 0     |

3. Example submission script below.

**Legend**:
**APPTAINER_CACHEDIR** = Sets the path to where `apptainer` will cache the containers used in the pipeline.
**input** = Path to `tsv` file which has the uniprot IDs and `bait` status.
**output** = Path to result directory.
**org_ref** = Path to the `tsv` file downloaded from Uniprot.
**num_recycles_colabfold** = Number of recycles to use in Colabfold. Higher the number better the prediction, but the slower the pipeline.
**resume** = Enables the pipeline to be used repeatedly. The pipeline will only run incomplete processes when rerun with the same inputs. 

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

module load nextflow/24.04.3 apptainer/latest

APPTAINER_CACHEDIR=/path/to/.apptainer/cache  # Path to cache directory for apptainer cache

export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

nextflow run baldikacti/chienlab-proteinfold \
      --input /path/to/acclist.tsv \
      --outdir /path/to/results \
      --org_ref /path/to/organism_reference.tsv \
      --mode colabfold \
      --num_recycles_colabfold 5 \
      --colabfold_model_preset "alphafold2_multimer_v3" \
      -profile unity \
      -resume
```

4. Submit to slurm with `sbatch main.sh`