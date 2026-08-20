/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH                       } from '../modules/colabfold_batch'
include { POOL                                  } from '../modules/pool'
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
        ch_input_raw = POOL.out.pool_fasta
            .flatten()
            .map { v -> tuple(v.getBaseName(), v) }
    } else {
        PROCESS_TSV (accession_file, 'colabfold')
        ch_input_raw = PROCESS_TSV.out.processed_tsv_output
            .flatten()
            .map { v -> tuple(v.getBaseName(), v) }
    }

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
    preprocessed = PROCESS_TSV.out.processed_tsv_output
    predictions  = COLABFOLD_BATCH.out.results
    ranked       = RANK_AF.out.tsv
}