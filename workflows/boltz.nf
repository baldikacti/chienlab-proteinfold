/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV           } from '../modules/process_tsv'
include { POOL                  } from '../modules/pool'
include { PREPARE_BOLTZ_CACHE   } from '../modules/prepare_boltz_cache'
include { BOLTZ_PREDICT         } from '../modules/boltz_predict'
include { RANK_AF               } from '../modules/rank_af'

workflow BOLTZ {
    take:
    accession_file
    boltz_model
    inf_batch

    main:

    if (params.pool) {
        // Pool mode input is a single-row, two-column TSV: <bait fasta>\t<pool fasta>
        ch_pool_input = accession_file
            .splitCsv(sep: '\t', strip: true)
            .first()
            .map { row ->
                if (row.size() < 2) {
                    error "Pool input '${params.input}' must have two tab-separated columns: bait FASTA and pool FASTA"
                }
                tuple(file(row[0], checkIfExists: true), file(row[1], checkIfExists: true))
            }

        POOL (ch_pool_input)
        ch_input_raw = POOL.out.pool_fasta.flatten()
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
    msa                 = BOLTZ_PREDICT.out.msa
    predictions         = BOLTZ_PREDICT.out.predictions
    ranked              = RANK_AF.out.tsv
}