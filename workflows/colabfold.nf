/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH                       } from '../modules/colabfold_batch'
include { COLABFOLD_BATCH as COLABFOLD_BATCH_TOP} from '../modules/colabfold_batch'
include { PREPARE_COLABFOLD_CACHE               } from '../modules/prepare_colabfold_cache'
include { PREPROCESS_FASTA_PAIRS                } from '../modules/preprocess_fasta_pairs'
include { RANK_IPTM                             } from '../modules/rank_iptm'


workflow COLABFOLD {
    take:
    accession_file
    colabfold_model_preset
    num_recycles_colabfold

    main:

    // Creates channel for the organism reference file if exists, otherwise uses a placeholder
    organism_ref = params.org_ref ? Channel.fromPath(params.org_ref, checkIfExists: true) : Channel.fromPath("$projectDir/assets/NO_FILE")

    //
    // Create input channel from input file provided through params.input
    //
    PREPROCESS_FASTA_PAIRS(accession_file)
    ch_fasta = PREPROCESS_FASTA_PAIRS.out.fasta
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

    RANK_IPTM(
        COLABFOLD_BATCH.out.json.collect(),
        organism_ref
        )
    ch_ranked = RANK_IPTM.out.tsv

    if (params.top_rank) {        

        ch_ranked_fasta = ch_ranked
            .splitCsv(header: true, sep: "\t", limit: params.top_rank)
            .map { tuple(it.paired_uniprotID) }

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