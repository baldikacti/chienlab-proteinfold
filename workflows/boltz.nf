/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV           } from '../modules/process_tsv'
include { POOL; poolInput       } from '../modules/pool'
include { PREPARE_BOLTZ_CACHE   } from '../modules/prepare_boltz_cache'
include { BOLTZ_PREDICT         } from '../modules/boltz_predict'
include { RANK_AF               } from '../modules/rank_af'

workflow BOLTZ {
    take:
    accession_file
    boltz_model
    inf_batch

    main:

    ch_pools_tsv = channel.empty()

    if (params.pool) {
        POOL (poolInput(accession_file), 'boltz')
        ch_input_raw = POOL.out.pools.flatten()
        ch_pools_tsv = POOL.out.pools_tsv
    } else {
        PROCESS_TSV (accession_file, 'boltz')
        ch_input_raw = PROCESS_TSV.out.processed_tsv_output.flatten()
    }

    PREPARE_BOLTZ_CACHE(boltz_model)
    boltz_cache = PREPARE_BOLTZ_CACHE.out.cache
    
    BOLTZ_PREDICT (
        ch_input_raw.collate( inf_batch ),
        boltz_cache
    )

    RANK_AF (
        BOLTZ_PREDICT.out.confidence_json.collect(),
        'boltz'
    )

    emit:
    preprocessed        = ch_input_raw
    pools_tsv           = ch_pools_tsv
    msa                 = BOLTZ_PREDICT.out.msa
    predictions         = BOLTZ_PREDICT.out.predictions
    ranked              = RANK_AF.out.tsv
}