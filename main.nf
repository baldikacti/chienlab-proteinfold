#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    baldikacti/chienlab-proteinfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/baldikacti/chienlab-proteinfold
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD             } from './workflows/colabfold'
include { ALPHAFOLD3            } from './workflows/alphafold3'
include { BOLTZ                 } from './workflows/boltz'
include { checkRequiredParams   } from './modules/checkparams.nf'


/*
================================================================================
    Workflow
================================================================================
*/


workflow {
    main:


    // Check and validate input paramaters
    checkRequiredParams()

    channel
        .fromPath(params.input, checkIfExists: true)
        .set { ch_input }

    ch_preprocessing        = channel.empty()
    ch_pools_tsv            = channel.empty()
    ch_af3_msa              = channel.empty()
    ch_af3_folds            = channel.empty()
    ch_boltz_folds          = channel.empty()
    ch_colabfold            = channel.empty()
    ch_ranked               = channel.empty()

    if (params.mode == "colabfold") {
        COLABFOLD (
            ch_input,
            params.num_recycle
        )

        ch_preprocessing        = COLABFOLD.out.preprocessed
        ch_pools_tsv            = COLABFOLD.out.pools_tsv
        ch_colabfold            = COLABFOLD.out.predictions
        ch_ranked               = COLABFOLD.out.ranked
    } else if (params.mode == "alphafold3") {

        ch_af3_db = file(params.db_dir, checkIfExists: true)
        ch_model_dir = file(params.model_dir, checkIfExists: true)

        ALPHAFOLD3 (
            ch_input,
            ch_af3_db,
            ch_model_dir,
            params.inf_batch
        )

        ch_preprocessing        = ALPHAFOLD3.out.preprocessed
        ch_pools_tsv            = ALPHAFOLD3.out.pools_tsv
        ch_af3_msa              = ALPHAFOLD3.out.msa
        ch_af3_folds            = ALPHAFOLD3.out.folds
        ch_ranked               = ALPHAFOLD3.out.ranked
    } else if (params.mode == "boltz") {

        BOLTZ (
            ch_input,
            params.model,
            params.inf_batch
        )

        ch_preprocessing        = BOLTZ.out.preprocessed
        ch_pools_tsv            = BOLTZ.out.pools_tsv
        ch_boltz_folds          = BOLTZ.out.msa.mix(BOLTZ.out.predictions)
        ch_ranked               = BOLTZ.out.ranked
    }

    publish:
    preprocessing           = ch_preprocessing
    pools_summary           = ch_pools_tsv
    alphafold3_msa          = ch_af3_msa
    alphafold3_folds        = ch_af3_folds
    boltz_folds             = ch_boltz_folds
    colabfold_predictions   = ch_colabfold
    ranked_results          = ch_ranked

}

/*
================================================================================
    Outputs
================================================================================
*/

// Paths are relative to `outputDir` (set from `params.outdir` in nextflow.config).
// The publish mode comes from `workflow.output.mode`.
output {
    preprocessing {
        path "${params.mode}/preprocessing"
    }

    // POOL emits `pools.tsv`, one row per generated pool
    pools_summary {
        path "${params.mode}/preprocessing"
    }

    // AF3_MSA emits `*_data.json`
    alphafold3_msa {
        path 'alphafold3/msa'
    }

    // AF3_FOLD emits `folds/*`
    alphafold3_folds {
        path 'alphafold3'
    }

    // BOLTZ_PREDICT emits `folds/msa/*` and `folds/predictions/*`
    boltz_folds {
        path 'boltz'
    }

    // COLABFOLD_BATCH emits every result file for one accession, published together
    colabfold_predictions {
        path { accID, files -> files >> "colabfold/${accID}/" }
    }

    ranked_results {
        path '.'
    }
}