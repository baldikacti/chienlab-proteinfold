process RANK_IPTM {
    label 'process_single'
    publishDir "${params.outdir}/ranked_results", mode: 'copy', pattern: '*.tsv'

    container "docker://baldikacti/chienlab_proteinfold_rverse:4.4.2"

    input:
    path json
    path org_ref

    output:
    path ("*.tsv") , emit: tsv

    script:
    def ref = org_ref.name != 'NO_FILE' ? "$org_ref" : 'NO_REF'
    """
    rank_pairs.R $ref
    """
}