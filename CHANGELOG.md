# Version v0.3.0

- Adds a new paramater, `top_rank`, that reruns the top selected number of predictions with 20 recycles and exports them to a new folder called `toprank`.
- Added CHANGELOG.md
- Cleans up the `nextflow.config` to use the `unity` profile on `nf-core/configs`.

# Version v0.2.0

- Fixes a bug that made each colabfold call to download the colabfold cache over and over again.

# Version v0.1.0

- Initial release
- Contains Colabfold pipeline for the Unity cluster
