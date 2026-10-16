# Pipeline preparation changes

## 16S.sh

- Removed explanatory comments and separator banners; retained shebang and literal # in data patterns.
- CLEAN_DIR is now required instead of using a private server directory.
- BASE_WORK defaults to results/qiime2 under the current directory; all original analysis settings remain unchanged.
- Retained existing protection that refuses an already-existing WORK directory.

## continue_ASV_taxonomy_tree.sh

- Removed explanatory comments and separator banners; retained shebang and literal # in export headers/patterns.
- WORK is required; SILVA_DIR, REF_SEQS, REF_TAX, CONDA_SH, CLASSIFIER_DIR and resource settings can be supplied through environment variables.
- Resolved relative input and output paths before changing directory.
- All generated files are written to a new RESULT_DIR (default WORK/taxonomy_tree), leaving the upstream run untouched.
- Removed every rm command. Existing RESULT_DIR is refused; invalid cached classifier triggers failure instead of deletion.
- Added early checks for the two frequency artifacts already required by later export steps.
- Removed unused THREADS variable. Classification and tree inference defaults are unchanged.

## run_icamp_bin48.R

- Removed explanatory comments and separator banners; retained the Rscript shebang.
- Replaced server paths with ASV_EXPORT_DIR, ASV_TABLE_FILE, TREE_FILE, ICAMP_R_LIBRARY and ICAMP_OUTPUT_DIR environment variables.
- Refuse a nonempty output directory to prevent mixing current results with an old phylogenetic distance cache.
- Scientific inputs, zero filtering, tree matching, seed and every icamp.big argument remain unchanged.

## mbpd_asv_direct.py

- Removed module docstring; retained the Python shebang; formatting normalized without changing BLAST arguments.
- Retained original default coverage 90.0. Manuscript reproduction example explicitly passes --coverage 95.
- Refuse a nonempty output directory to prevent reuse of a stale BLAST database or overwrite of previous results.
- Fixed empty FASTA identifiers to raise the intended ValueError instead of an uncaught IndexError.
