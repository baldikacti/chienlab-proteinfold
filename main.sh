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
ACC_FILE=$WORKDIR/tests/acclist.csv                                   # Path to bait:prey csv file
CCNAREF=$WORKDIR/references/CCNA_ref.csv                              # Path to gene reference csv file
UNIPROT_MAP=$WORKDIR/references/cc_uniprot_mappings.csv               # Path to UniprotID:GeneID mapping csv file
APPTAINER_CACHEDIR=/work/pi_pchien_umass_edu/berent/.apptainer/cache  # Path to cache directory for apptainer cache

cd $WORKDIR || exit

# DO NOT CHANGE THESE VARIABLES
SAMPLE_SHEET=$RESULTDIR/samplesheet.csv
FASTA_DIR=$RESULTDIR/fasta
WORK_BIN=$WORKDIR/bin
R_CONTAINER=$WORK_BIN/chienlab_colabfold_4.4.sif
COLABFOLD_CACHE=$WORKDIR/cache

export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

## Download fasta files from uniprotIDs and generate bait:pair fasta files

if [ ! -f "$R_CONTAINER" ]; then
  echo "Pulling container baldikacti/chienlab_colabfold:4.4"
  apptainer pull --dir "$WORK_BIN" docker://baldikacti/chienlab_colabfold:4.4
fi

if [ ! -d "$FASTA_DIR" ] || [ ! -f "$SAMPLE_SHEET" ]; then
  $R_CONTAINER Rscript "$WORK_BIN"/combine_fasta.R \
    --acc_file "$ACC_FILE" \
    --fasta_dir "$FASTA_DIR" \
    --samplesheet "$SAMPLE_SHEET"
fi

if [ ! -d "$COLABFOLD_CACHE" ]; then
    mkdir -p $COLABFOLD_CACHE && cd $COLABFOLD_CACHE
    wget https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar
    tar -xavf alphafold_params_colab_2022-12-06.tar
    rm alphafold_params_colab_2022-12-06.tar
    cd $WORKDIR
fi


## Run ColabFold using nf-core/proteinfold pipeline ####

nextflow run main.nf \
      --input "$SAMPLE_SHEET" \
      --outdir "$RESULTDIR" \
      --mode colabfold \
      --colabfold_cache "$COLABFOLD_CACHE" \
      --num_recycles_colabfold 3 \
      --colabfold_model_preset "alphafold2_multimer_v3" \
      -c conf/unity.config \
      -profile unity \
      -resume

# Sumamarise and export ranked results to csv file
$R_CONTAINER Rscript "$WORK_BIN"/rank_pairs.R "$RESULTDIR" "$CCNAREF" "$UNIPROT_MAP"