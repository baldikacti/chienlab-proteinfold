process PREPARE_BOLTZ_CACHE {
    label "process_single"

    container "docker://baldikacti/boltz:latest"

    input:
    val model

    output:
    path ("cache/")         , emit: cache

    script:
    """
    prepare_boltz_cache.py ./cache $model
    """
}