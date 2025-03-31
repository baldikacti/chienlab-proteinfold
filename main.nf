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
    Validate inputs and create channels for files
================================================================================
*/

if (params.input) {
    ch_input = file(params.input, checkIfExists: true)
} else { exit 1, 'Input file not specified!' }

if (params.org_ref) {
    ch_org_ref = file(params.org_ref, checkIfExists: true)
} else { exit 1, 'Organism reference file not specified!' }


/*
================================================================================
    Workflow
================================================================================
*/

workflow {
    main:


    COLABFOLD (
        ch_input,
        params.colabfold_model_preset,
        params.num_recycles_colabfold,
        ch_org_ref
    )

}