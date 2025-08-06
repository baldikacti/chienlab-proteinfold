/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH                       } from '../modules/colabfold_batch'
include { COLABFOLD_BATCH as COLABFOLD_BATCH_TOP} from '../modules/colabfold_batch'
include { PREPARE_COLABFOLD_CACHE               } from '../modules/prepare_colabfold_cache'
include { PROCESS_TSV                           } from '../modules/process_tsv'
include { RANK_AF                               } from '../modules/rank_af'


workflow COLABFOLD {
    take:
    accession_file
    colabfold_model_preset
    num_recycles_colabfold

    main:

    //
    // Create input channel from input file provided through params.input
    //
    PROCESS_TSV(accession_file, 'colabfold')
    ch_fasta = PROCESS_TSV.out.processed_tsv_output
        .flatten()
        .map { tuple(it.getBaseName(), it) }

    PREPARE_COLABFOLD_CACHE()
    colabfold_cache = PREPARE_COLABFOLD_CACHE.out.cache

    COLABFOLD_BATCH(
            ch_fasta,
            colabfold_model_preset,
            colabfold_cache,
            num_recycles_colabfold,
            "screen"
        )

    RANK_AF(
        COLABFOLD_BATCH.out.json.collect(),
        'colabfold'
        )
    ch_ranked = RANK_AF.out.tsv

    if (params.top_rank) {        

        ch_ranked_fasta = ch_ranked
            .splitCsv(header: true, sep: "\t", limit: params.top_rank)
            .map { tuple(it.foldid) }

        // Filter ch_fasta based on ch_ranked_fasta
        ch_filtered_fasta = ch_fasta
            .join(ch_ranked_fasta)
        
        COLABFOLD_BATCH_TOP(
            ch_filtered_fasta,
            colabfold_model_preset,
            colabfold_cache,
            20,
            "toprank"
        )
    }
}