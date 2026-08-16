process AFPOOL {
    label "process_single"

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    path bait_fasta
    path pool_fasta

    output:
    path ("results/pool_*.fasta")         , emit: pool_fasta

    script:
    def args = task.ext.args ?: ''
    """
    afpool.py --bait ${bait_fasta} --pool ${pool_fasta} --output ./results $args
    """
}