/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH        } from '../modules/colabfold_batch'
include { PREPARE_COLABFOLD_CACHE} from '../modules/prepare_colabfold_cache'
include { PREPROCESS_FASTA_PAIRS } from '../modules/preprocess_fasta_pairs'


workflow COLABFOLD {
    take:
    accession_file
    colabfold_model_preset
    num_recycles_colabfold

    main:
    //
    // Create input channel from input file provided through params.input
    //
    PREPROCESS_FASTA_PAIRS(accession_file)
    ch_fasta = PREPROCESS_FASTA_PAIRS.out.fasta.flatten()

    //Channel.fromPath(PREPROCESS_FASTA_PAIRS.out.processed_input)
    //    .splitCsv(header: true)
    //    .map { row-> tuple(row.sequence, file(row.fasta)) }
    //    .set { ch_fasta }

    PREPARE_COLABFOLD_CACHE()

    COLABFOLD_BATCH(
            ch_fasta,
            colabfold_model_preset,
            PREPARE_COLABFOLD_CACHE.out.cache,
            num_recycles_colabfold
        )
}