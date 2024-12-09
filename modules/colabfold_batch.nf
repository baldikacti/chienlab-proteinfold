process COLABFOLD_BATCH {
    tag "$meta"
    label 'process_medium'
    label "gpu"
    publishDir = [
                path: { "${params.outdir}/${params.mode}/${meta}" },
                mode: 'copy',
                pattern: '*.*'
            ]

    container "docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"

    input:
    tuple val(meta), path(fasta)
    val cb_model
    path cb_cache
    val numRec

    output:
    path ("*")         , emit: pdb
    path ("*.png") , emit: multiqc

    script:
    """
    colabfold_batch \\
        ${fasta} \\
        \$PWD \\
        --num-recycle ${numRec} \\
        --msa-mode 'mmseqs2_uniref_env' \\
        --model-type ${cb_model} \\
        --templates \\
        --amber \\
        --use-gpu-relax \\
        --data ${cb_cache} \\
        --host-url http://cfold-db:8888
    """
}