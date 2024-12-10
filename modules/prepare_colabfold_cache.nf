process PREPARE_COLABFOLD_CACHE {
    label "process_single"

    container "docker://baldikacti/chienlab_proteinfold_ubuntu:noble"

    output:
    path ("*")         , emit: cache

    script:
    """
    wget --no-check-certificate https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar
    tar -xavf alphafold_params_colab_2022-12-06.tar
    rm alphafold_params_colab_2022-12-06.tar
    """
}