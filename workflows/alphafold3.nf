/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS_AF3    } from '../modules/preprocess_af3'
include { AF3_MSA           } from '../modules/af3_msa'
include { AF3_FOLD          } from '../modules/af3_fold'
include { AF3_RANK          } from '../modules/af3_rank'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow ALPHAFOLD3 {
    take:
    accession_file
    database_dir
    model_dir

    main:

    PREPROCESS_AF3 (
        accession_file
    )
    ch_json_raw = PREPROCESS_AF3.out.json.flatten()

    AF3_MSA (
        ch_json_raw,
        database_dir
    )
    msa_json = AF3_MSA.out.af3_json_processed

    ch_msa_json = msa_json
        .collate( params.inf_batch )

    AF3_FOLD (
        ch_msa_json,
        database_dir,
        model_dir
    )
    ch_json_confidence = AF3_FOLD.out.summary_json.collect()

    AF3_RANK (
        ch_json_confidence
    )
}
