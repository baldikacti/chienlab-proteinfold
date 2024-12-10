process RANK_IPTM {
    label 'process_single'
    publishDir = [
                path: { "${params.outdir}/ranked_results" },
                mode: 'copy',
                pattern: '*.tsv'
            ]

    container "docker://baldikacti/chienlab_colabfold:4.4"

    input:
    path json

    output:
    path ("*.tsv") , emit: tsv

    script:
    """
    rank_pairs.R ${json}
    """
}