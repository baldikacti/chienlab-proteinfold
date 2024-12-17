#!/usr/bin/bash
#SBATCH --job-name=chienlab-proteinfold     # Job name
#SBATCH --partition=cpu                     # Partition (queue) name
#SBATCH -c 2                                # Number of CPUs
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --mem=10gb                          # Job memory request
#SBATCH --time=1-00:00:00                  # Time limit days-hrs:min:sec
##SBATCH -q long
#SBATCH --output=logs/chienlab-proteinfold_%j.log

module load nextflow/24.04.3 apptainer/latest

APPTAINER_CACHEDIR=/work/pi_pchien_umass_edu/berent/.apptainer/cache  # Path to cache directory for apptainer cache

export APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR
export NXF_OPTS="-Xms1G -Xmx8G"

nextflow run baldikacti/chienlab-proteinfold -r v0.2.0 \
      --input /work/pi_pchien_umass_edu/berent/chienlab-proteinfold/tests/acclist.tsv \
      --outdir /work/pi_pchien_umass_edu/berent/chienlab-proteinfold/results \
      --org_ref /work/pi_pchien_umass_edu/berent/chienlab-proteinfold/tests/uniprotkb_proteome_UP000001364_cc.tsv \
      --mode colabfold \
      --num_recycles_colabfold 5 \
      --colabfold_model_preset "alphafold2_multimer_v3" \
      -profile unity \
      -resume