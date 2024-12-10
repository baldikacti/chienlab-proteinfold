#!/usr/bin/bash
#SBATCH --job-name=chienlab-proteinfold       # Job name
#SBATCH --partition=cpu                       # Partition (queue) name
#SBATCH -c 2                                  # Number of CPUs
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --mem=34gb                            # Job memory request
#SBATCH --time=2-00:00:00                    # Time limit hrs:min:sec
##SBATCH -q long
#SBATCH --output=logs/chienlab-proteinfold_%j.log

set -eou pipefail

module load nextflow/24.04.3 apptainer/latest

# CHANGE THE PATHS BELOW
WORKDIR=/work/pi_pchien_umass_edu/berent/chienlab-proteinfold         # Path to pipeline location (Should be under a fast filesystem like '/work' or '/scratch')
RESULTDIR=$WORKDIR/results/tests                                      # Path to results
ACC_FILE=$WORKDIR/tests/acclist.tsv                                   # Path to bait:prey csv file
APPTAINER_CACHEDIR=/work/pi_pchien_umass_edu/berent/.apptainer/cache  # Path to cache directory for apptainer cache

cd $WORKDIR || exit

# DO NOT CHANGE THESE VARIABLES
export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

## Run ColabFold using nf-core/proteinfold pipeline ####

nextflow run main.nf \
      --input "$ACC_FILE" \
      --outdir "$RESULTDIR" \
      --mode colabfold \
      --num_recycles_colabfold 3 \
      --colabfold_model_preset "alphafold2_multimer_v3" \
      -profile unity,debug \
      -resume
