process AF3_RANK {
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy', pattern: 'ranked_results.tsv'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    path summary_json

    output:
    path ("ranked_results.tsv")

    script:
    """
    rank_af3.py --output=ranked_results.tsv
    """
}