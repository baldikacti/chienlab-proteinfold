process BOLTZ_PREDICT {
    label 'gpu'
    label 'error_ignore'
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "folds/**"

    container "docker://baldikacti/boltz:latest"

    input:
    path ("input_fasta/*")
    path cache

    output:
    path ("folds/**")
    path ("folds/predictions/*/*_model_0.json"), emit: confidence_json

    script:
    def args = task.ext.args ?: ''
    """
    mkdir numba_cache
    mkdir pytorch_kernel_cache
    export NUMBA_CACHE_DIR=./numba_cache
    export PYTORCH_KERNEL_CACHE_PATH=./pytorch_kernel_cache
    boltz predict --out_dir . --cache $cache --num_workers $task.cpus $args input_fasta/
    mv boltz_results_input_fasta folds
    """
}
