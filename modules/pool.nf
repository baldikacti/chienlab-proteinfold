process POOL {
    label 'process_single'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    tuple path(bait_file), path(pool_file)

    output:
    path ("pools/*.fasta") , emit: pool_fasta

    script:
    """
    afpool.py -b $bait_file -p $pool_file -o pools
    """
}