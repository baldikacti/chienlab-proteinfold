process PREPROCESS_FASTA_PAIRS {
    label 'process_single'
    publishDir = [
                path: { "${params.outdir}/fasta" },
                mode: 'copy',
                pattern: '*.fasta'
            ]

    container "docker://baldikacti/chienlab_colabfold:4.4"

    input:
    path acc_file

    output:
    path ("*.fasta") , emit: fasta

    script:
    """
    combine_fasta.R --acc_file ${acc_file}
    """
}