process PREPROCESS_FASTA_PAIRS {
    label 'process_single'
    publishDir = [
                path: { "${params.outdir}/fasta" },
                mode: 'copy',
                pattern: '*.fasta'
            ]

    container "docker://baldikacti/chienlab_proteinfold_rverse:4.4.2"

    input:
    path acc_file

    output:
    path ("*.fasta") , emit: fasta

    script:
    """
    combine_fasta.R --acc_file ${acc_file}
    """
}