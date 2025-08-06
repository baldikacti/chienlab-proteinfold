process RANK_AF {
    label 'process_single'
    publishDir "${params.outdir}", mode: 'copy', pattern: "*ranked_results.tsv"

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    path summary_json
    val mode

    output:
    path ("*ranked_results.tsv"), emit: tsv

    script:
    """
    rank_af.py --output="${mode}_ranked_results.tsv" --mode $mode
    """
}