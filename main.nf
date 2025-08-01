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
include { checkRequiredParams   } from './lib/checkparams.groovy'


/*
================================================================================
    Workflow
================================================================================
*/


workflow {
    main:


    // Check and validate input paramaters
    checkRequiredParams()

    Channel
        .fromPath(params.input, checkIfExists: true)
        .set { ch_input }

    if (params.mode == "colabfold") {
        COLABFOLD (
            ch_input,
            params.colabfold_model_preset,
            params.num_recycles_colabfold
        )
    } else if (params.mode == "alphafold3") {
        
        ch_af3_db = file(params.db_dir, checkIfExists: true)
        ch_model_dir = file(params.model_dir, checkIfExists: true)
        
        ALPHAFOLD3 (
            ch_input,
            ch_af3_db,
            ch_model_dir
        )
    }

}