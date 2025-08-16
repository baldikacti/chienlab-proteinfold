process COLABFOLD_BATCH {
    tag "$accID"
    label 'gpu'
    label 'error_ignore'
    publishDir "${params.outdir}/${params.mode}/$outDir/$accID", mode: 'copy', pattern: '*.*'

    container "docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"

    input:
    tuple val(accID), path(fasta)
    path ("params/*")
    val  numRec
    val  outDir

    output:
    path ("*")                        , emit: pdb
    path ("*_toprank.json")           , emit: json
    path ("*.png")                    , emit: multiqc

    script:
    def args = task.ext.args ?: ''
    """
    colabfold_batch \\
        ${fasta} \\
        \$PWD \\
        --num-recycle ${numRec} \\
        --data \$PWD \\
        $args
    
    ln -s *scores_rank_001_*.json ${accID}_toprank.json
    """
}