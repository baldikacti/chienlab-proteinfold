/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_TSV       } from '../modules/process_tsv'
include { POOL              } from '../modules/pool'
include { FASTA2JSON        } from '../modules/fasta2json'
include { AF3_MSA           } from '../modules/af3_msa'
include { AF3_FOLD          } from '../modules/af3_fold'
include { RANK_AF           } from '../modules/rank_af'

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
        ch_pool_raw = POOL.out.pool_fasta.flatten()

        FASTA2JSON (ch_pool_raw)
        ch_input_raw = FASTA2JSON.out.af3_json.flatten()
    } else {
        PROCESS_TSV (accession_file, 'alphafold3')
        ch_input_raw = PROCESS_TSV.out.processed_tsv_output.flatten()
    }

    AF3_MSA (
        ch_input_raw,
        database_dir
    )
    msa_json = AF3_MSA.out.af3_json_processed

    ch_msa_json = msa_json
        .collate( inf_batch )

    AF3_FOLD (
        ch_msa_json,
        database_dir,
        model_dir
    )
    ch_json_confidence = AF3_FOLD.out.summary_json.collect()

    RANK_AF (
        ch_json_confidence,
        'alphafold3'
    )

    emit:
    preprocessed        = ch_input_raw
    preprocessed_pool   = ch_pool_raw
    msa                 = AF3_MSA.out.af3_json_processed
    folds               = AF3_FOLD.out.folds
    ranked              = RANK_AF.out.tsv
}
