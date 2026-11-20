/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH                       } from '../modules/colabfold_batch'
include { POOL; poolInput                       } from '../modules/pool'
include { PREPARE_COLABFOLD_CACHE               } from '../modules/prepare_colabfold_cache'
include { PROCESS_TSV                           } from '../modules/process_tsv'
include { RANK_AF                               } from '../modules/rank_af'


workflow COLABFOLD {
    take:
    accession_file
    num_recycle

    main:

    //
    // Create input channel from input file provided through params.input
    //
    ch_pools_tsv = channel.empty()

    if (params.pool) {
        POOL (poolInput(accession_file), 'colabfold')
        ch_preprocessed = POOL.out.pools.flatten()
        ch_pools_tsv = POOL.out.pools_tsv
    } else {
        PROCESS_TSV (accession_file, 'colabfold')
        ch_preprocessed = PROCESS_TSV.out.processed_tsv_output.flatten()
    }

    ch_input_raw = ch_preprocessed.map { v -> tuple(v.getBaseName(), v) }

    PREPARE_COLABFOLD_CACHE()
    colabfold_cache = PREPARE_COLABFOLD_CACHE.out.cache

    COLABFOLD_BATCH(
            ch_input_raw,
            colabfold_cache,
            num_recycle
        )

    RANK_AF(
        COLABFOLD_BATCH.out.json.collect(),
        'colabfold'
        )

    emit:
    preprocessed = ch_preprocessed
    pools_tsv    = ch_pools_tsv
    predictions  = COLABFOLD_BATCH.out.results
    ranked       = RANK_AF.out.tsv
}