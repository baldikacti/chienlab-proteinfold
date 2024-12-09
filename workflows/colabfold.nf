/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COLABFOLD_BATCH        } from '../modules/colabfold_batch'
//include { fromSamplesheet        } from 'plugin/nf-validation'


workflow COLABFOLD {
    take:
    colabfold_model_preset
    colabfold_cache
    num_recycles_colabfold

    main:
    //
    // Create input channel from input file provided through params.input
    //
    Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row-> tuple(row.sequence, file(row.fasta)) }
        .set { ch_fasta }
    //Channel
    //    .fromSamplesheet("input")
    //    .set { ch_fasta }

    COLABFOLD_BATCH(
            ch_fasta,
            colabfold_model_preset,
            colabfold_cache,
            num_recycles_colabfold
        )
}