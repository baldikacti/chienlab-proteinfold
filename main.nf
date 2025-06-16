#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    baldikacti/chienlab-proteinfold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/baldikacti/chienlab-proteinfold
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD             } from './workflows/colabfold'


/*
================================================================================
    Workflow
================================================================================
*/


workflow {
    main:


    /*
        Validate inputs and create channels for files
    */
    if (!params.input) {
        exit 1, 'Input file not specified!'
    }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .set { ch_input }

    // Creates channel for the organism reference file if exists, otherwise uses a placeholder
    ch_org_ref = params.org_ref ? Channel.fromPath(params.org_ref, checkIfExists: true) : Channel.fromPath("$projectDir/assets/NO_FILE")


    COLABFOLD (
        ch_input,
        params.colabfold_model_preset,
        params.num_recycles_colabfold,
        ch_org_ref
    )
}