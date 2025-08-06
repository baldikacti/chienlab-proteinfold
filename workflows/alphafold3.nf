/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV       } from '../modules/process_tsv'
include { AF3_MSA           } from '../modules/af3_msa'
include { AF3_FOLD          } from '../modules/af3_fold'
include { RANK_AF           } from '../modules/rank_af'

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

    PROCESS_TSV (accession_file, 'alphafold3')
    ch_json_raw = PROCESS_TSV.out.processed_tsv_output.flatten()

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

    RANK_AF (
        ch_json_confidence,
        'alphafold3'
    )
}
