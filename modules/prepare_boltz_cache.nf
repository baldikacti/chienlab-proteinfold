process PREPARE_BOLTZ_CACHE {
    label "process_single"

    container "docker://registry.harmonyhpc.io/baldikacti_umass_edu-pub/chienlab_proteinfold/boltz:latest"

    input:
    val model

    output:
    path ("cache/")         , emit: cache

    script:
    """
    prepare_boltz_cache.py ./cache $model
    """
}