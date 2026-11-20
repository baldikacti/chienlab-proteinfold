/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV       } from '../modules/process_tsv'
include { POOL; POOL_ANALYSIS; poolInput   } from '../modules/pool'
include { AF3_MSA           } from '../modules/af3_msa'
include { AF3_FOLD          } from '../modules/af3_fold'
include { RANK_AF           } from '../modules/rank_af'
include { batchByInput      } from '../modules/utils'

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
    inf_batch

    main:

    ch_pools_tsv = channel.empty()

    if (params.pool) {
        POOL (poolInput(accession_file), 'alphafold3')
        ch_input_raw = POOL.out.pools.flatten()
        ch_pools_tsv = POOL.out.pools_tsv
    } else {
        PROCESS_TSV (accession_file, 'alphafold3')
        ch_input_raw = PROCESS_TSV.out.processed_tsv_output.flatten()
    }

    AF3_MSA (
        ch_input_raw,
        database_dir
    )
    msa_json = AF3_MSA.out.af3_json_processed

    AF3_FOLD (
        batchByInput(ch_input_raw, msa_json, inf_batch, '_data'),
        database_dir,
        model_dir
    )
    ch_json_confidence = AF3_FOLD.out.summary_json.collect()
    ch_af3_dir = AF3_FOLD.out.folds

    if (params.pool) {
        POOL_ANALYSIS (
            ch_af3_dir.collect(),
            ch_pools_tsv
        )
        af3_results = POOL_ANALYSIS.out.analysis_results
    } else {
        RANK_AF (
            ch_json_confidence,
            'alphafold3'
        )
        af3_results = RANK_AF.out.tsv
    }

    emit:
    preprocessed        = ch_input_raw
    pools_tsv           = ch_pools_tsv
    msa                 = AF3_MSA.out.af3_json_processed
    folds               = AF3_FOLD.out.folds
    ranked              = af3_results
}
