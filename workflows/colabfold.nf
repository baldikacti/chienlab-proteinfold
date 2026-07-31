/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH                       } from '../modules/colabfold_batch'
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
    PROCESS_TSV(accession_file, 'colabfold')
    ch_fasta = PROCESS_TSV.out.processed_tsv_output
        .flatten()
        .map { v -> tuple(v.getBaseName(), v) }

    PREPARE_COLABFOLD_CACHE()
    colabfold_cache = PREPARE_COLABFOLD_CACHE.out.cache

    COLABFOLD_BATCH(
            ch_fasta,
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