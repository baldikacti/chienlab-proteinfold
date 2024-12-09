process PREPARE_COLABFOLD_CACHE {
    label "process_single"

    container "ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"

    output:
    path ("cache/*")         , emit: cache

    script:
    """
    mkdir ./cache
    python -m colabfold.download
    """
}