process COLABFOLD_BATCH {
    tag "$fastaID"
    label 'process_medium'
    label "gpu"
    publishDir = [
                path: { "${params.outdir}/${params.mode}/$fastaID" },
                mode: 'copy',
                pattern: '*.*'
            ]

    container "docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"

    input:
    path fasta
    val  cb_model
    path ("params/*")
    val  numRec

    output:
    path ("*")                  , emit: pdb
    path ("*scores_rank*.json") , emit: json
    path ("*.png")              , emit: multiqc

    script:
    fastaID = fasta.getSimpleName()
    
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
        --data \$PWD \\
        --host-url http://cfold-db:8888
    """
}