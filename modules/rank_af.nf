process RANK_AF {
    label 'process_single'

    container "docker://registry.harmonyhpc.io/baldikacti_umass_edu-pub/chienlab_proteinfold/chienlab_proteinfold_py:latest"

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