process RANK_IPTM {
    label 'process_single'
    publishDir = [
                path: { "${params.outdir}/ranked_results" },
                mode: 'copy',
                pattern: '*.tsv'
            ]

    container "docker://baldikacti/chienlab_proteinfold_rverse:4.4.2"

    input:
    path json
    path org_ref

    output:
    path ("*.tsv") , emit: tsv

    script:
    """
    rank_pairs.R ${org_ref}
    """
}