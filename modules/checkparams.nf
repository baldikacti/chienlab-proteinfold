// Check required parameters are provided
def checkRequiredParams() {

    def missing = []
    if (params.mode == "colabfold") {
        if (!params.input) missing.add('--input')
        if (!params.outdir) missing.add('--outdir')
    } else if (params.mode == "alphafold3") {
        if (!params.input) missing.add('--input')
        if (!params.model_dir) missing.add('--model_dir') 
        if (!params.outdir) missing.add('--outdir')
        if (!params.db_dir) missing.add('--db_dir')
    } else if (params.mode == "boltz") {
        if (!params.input) missing.add('--input')
        if (!params.outdir) missing.add('--outdir')
        if (!params.model) missing.add('--model')
    } else {
        error "Either missing or incorrect paramater passed to `--mode`. Options: `alphafold3`, `colabfold`, or 'boltz'."
    }

    if (params.pool) {
        if (!params.pool_mode) {
            missing.add('--pool_mode')
        } else if (!(params.pool_mode in ['bait_vs_all', 'all_vs_all'])) {
            error "Incorrect paramater passed to `--pool_mode`. Options: `bait_vs_all` or `all_vs_all`."
        }
        if (params.init_pool && params.pool_mode != 'all_vs_all') {
            error "`--init_pool` is only used in `all_vs_all` pool mode."
        }
        if (params.max_pools && params.pool_mode != 'all_vs_all') {
            error "`--max_pools` is only used in `all_vs_all` pool mode."
        }
    }

    if (missing.size() > 0) {
        log.error "Missing required parameters: ${missing.join(', ')}"
        exit 1
    }
}