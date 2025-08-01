// Check required parameters are provided
def checkRequiredParams() {

    def missing = []
    if (params.mode == "colabfold") {
        if (!params.input) missing.add('--input')
        if (!params.colabfold_model_preset) missing.add('--colabfold_model_preset') 
        if (!params.outdir) missing.add('--outdir')
    } else if (params.mode == "alphafold3") {
        if (!params.input) missing.add('--input')
        if (!params.model_dir) missing.add('--model_dir') 
        if (!params.outdir) missing.add('--outdir')
        if (!params.db_dir) missing.add('--db_dir')
    } else {
        error "Either missing or incorrect paramater passed to `--mode`. Options: `alphafold3` or `colabfold`."
    }

    if (missing.size() > 0) {
        log.error "Missing required parameters: ${missing.join(', ')}"
        exit 1
    }
}