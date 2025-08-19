/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV           } from '../modules/process_tsv'
include { PREPARE_BOLTZ_CACHE   } from '../modules/prepare_boltz_cache'
include { BOLTZ_PREDICT         } from '../modules/boltz_predict'
include { RANK_AF               } from '../modules/rank_af'

workflow BOLTZ {
    take:
    ch_input
    boltz_model

    main:

    PROCESS_TSV(ch_input, 'boltz')
    ch_fasta = PROCESS_TSV.out.processed_tsv_output.flatten()

    PREPARE_BOLTZ_CACHE(boltz_model)
    boltz_cache = PREPARE_BOLTZ_CACHE.out.cache
    
    BOLTZ_PREDICT (
        ch_fasta.collate( params.inf_batch ),
        boltz_cache
    )

    RANK_AF (
        BOLTZ_PREDICT.out.confidence_json.collect(),
        'boltz'
    )
}