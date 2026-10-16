# Danjiangkou 16S analysis

All scripts are in `analysis/`. Numbered folders follow the analysis order; calculation and plotting scripts for the same analysis are stored together. `00_runtime` contains the launcher, workflow table and dependencies.

## Analysis order

| Folder | Analysis |
| --- | --- |
| [01_sequence_processing](analysis/01_sequence_processing/) | Paired-end reads, ASV inference, taxonomy and rooted tree |
| [02_taxonomy_screening](analysis/02_taxonomy_screening/) | Taxonomy tables and MBPD candidate screening |
| [03_data_preparation](analysis/03_data_preparation/) | Whole-community rarefaction and candidate extraction |
| [04_hydrology_environment](analysis/04_hydrology_environment/) | Discharge preparation, SRI/SDFAI and environment curves |
| [05_study_area](analysis/05_study_area/) | GIS input checks and sampling map |
| [06_community_composition](analysis/06_community_composition/) | Venn diagrams, taxonomy network and composition |
| [07_alpha_diversity](analysis/07_alpha_diversity/) | Alpha diversity, group comparisons and richness trends |
| [08_beta_diversity](analysis/08_beta_diversity/) | PCoA, group comparisons and temporal trends |
| [09_differential_abundance](analysis/09_differential_abundance/) | LEfSe |
| [10_distance_decay](analysis/10_distance_decay/) | Environmental, spatial, temporal and water-level decay |
| [11_environmental_ordination](analysis/11_environmental_ordination/) | RDA calculation and plots |
| [12_mantel_correlation](analysis/12_mantel_correlation/) | Mantel tests and environmental correlations |
| [13_variation_partitioning](analysis/13_variation_partitioning/) | Environment, time and space partitioning |
| [14_random_forest](analysis/14_random_forest/) | 19-predictor random forest excluding Q |
| [15_water_level_regression](analysis/15_water_level_regression/) | Community and candidate richness regressions |
| [16_community_assembly](analysis/16_community_assembly/) | iCAMP and assembly-process trends |
| [17_static_networks](analysis/17_static_networks/) | Static networks, topology and node removal |
| [18_temporal_networks](analysis/18_temporal_networks/) | Sliding-window networks and temporal trends |
| [19_cyanobacteria_dynamics](analysis/19_cyanobacteria_dynamics/) | Cyanobacterial extraction and monthly abundance |
| [20_cross_correlation](analysis/20_cross_correlation/) | CCF calculations, plots and panel assembly |
| [21_candidate_dynamics](analysis/21_candidate_dynamics/) | Candidate abundance and richness trends |
| [22_indicator_thresholds](analysis/22_indicator_thresholds/) | TITAN2, threshold curves, ridgelines and aligned panels |
| [23_result_exports](analysis/23_result_exports/) | Composition and diversity result tables |

Within a folder, run calculation scripts before plots. Files marked `support` in the workflow table provide functions, configuration or dependent exports. Files marked `alternative` are optional entry points. `14_random_forest/00_reference_random_forest.R` is a reference dependency; use `14.01` for the random-forest calculation.

## Requirements

- R 4.5.1; package versions are listed in [R_packages.csv](analysis/00_runtime/R_packages.csv).
- Python 3.9 or later; [requirements.txt](analysis/00_runtime/requirements.txt) lists packages used for workbook extraction and PDF layout.
- Linux, Bash and QIIME 2 2026.4 for sequence processing.
- NCBI BLAST+ (`blastn` and `makeblastdb`) for MBPD screening.

## Run an analysis

Run commands from the repository root. Keep the study input dataset in a separate working directory, preserving its existing filenames, columns and relative paths. The launcher selects each analysis's data working directory from [workflow.tsv](analysis/00_runtime/workflow.tsv). Scripts write generated results into the data working directory.

```bash
python analysis/00_runtime/01_run_analysis.py --list
python analysis/00_runtime/01_run_analysis.py --list --step 07
python analysis/00_runtime/01_run_analysis.py 03.01 --data-root /path/to/study_data --check
python analysis/00_runtime/01_run_analysis.py 03.01 --data-root /path/to/study_data
```

`--check` checks script and working-directory availability without running the analysis. Use `--rscript /path/to/Rscript` if Rscript is not on PATH. Pass script-specific arguments after `--`. `--working-directory` overrides the selected data workspace.

The launcher sets `DJK_CODE_ROOT` to this repository and `DJK_PROJECT_ROOT` to the data directory. To reuse the same data root for several commands, set `DJK_PROJECT_ROOT` once:

```bash
export DJK_PROJECT_ROOT=/path/to/study_data
```

In PowerShell, use `$env:DJK_PROJECT_ROOT = 'C:/path/to/study_data'`.

ASV tables have ASVs in rows and samples in columns. Sample identifiers must match the environment and group tables. Use original counts for shared-ASV, radial-taxonomy and iCAMP analyses. Diversity analyses use the complete-community table rarefied once to 21,000 reads per sample; candidate richness uses a subset of that same table. Input data and reference databases are supplied separately.

## Sequence processing

The main workflow expects `sample_R1.fastq.gz` and `sample_R2.fastq.gz` files. Set `ENV_NAME` and `CONDA_SH` if the QIIME environment or Conda installation differs from the script defaults. Use new output directories.

```bash
CLEAN_DIR=/path/to/primer_trimmed_reads \
BASE_WORK=/path/to/results/qiime2 \
RUN_TAG=analysis \
CLASSIFIER_QZA=/path/to/silva_classifier.qza \
bash analysis/01_sequence_processing/01_process_amplicons.sh
```

The second script resumes taxonomy and tree construction from existing `02_dada2` and `03_qc` artifacts. It is optional when the main workflow already completed these steps.

```bash
WORK=/path/to/results/qiime2/ASV_QIIME2_analysis \
REF_SEQS=/path/to/silva_sequences.qza \
REF_TAX=/path/to/silva_taxonomy.qza \
RESULT_DIR=/path/to/results/taxonomy_tree \
bash analysis/01_sequence_processing/02_assign_taxonomy_build_tree.sh

python analysis/02_taxonomy_screening/02_screen_candidate_asvs.py \
  --asv-fasta /path/to/rep_seqs/dna-sequences.fasta \
  --db-fasta /path/to/MBPD_reference.fasta \
  --db-tax /path/to/MBPD_taxonomy.tsv \
  --outdir /path/to/results/mbpd \
  --identity 97 --coverage 95 --threads 8 --max-target-seqs 50

ASV_EXPORT_DIR=/path/to/06_export \
ICAMP_OUTPUT_DIR=/path/to/results/icamp \
Rscript analysis/16_community_assembly/01_run_icamp.R
```

MBPD requires `--coverage 95` for the 95% filter; its default is 90. Its E-value threshold is `1e-20`. Optional `--abundance` accepts a comma-delimited ASV table. iCAMP uses original counts, seed 123, 1,000 randomizations, minimum bin size 48, 16 workers and a 50-GB memory allowance. `ASV_TABLE_FILE` and `TREE_FILE` override its individual input files.

## Additional arguments

These commands assume `DJK_PROJECT_ROOT` is set. Hydrological preparation runs without a script argument, followed by the monthly grid, Gamma-MLE SRI/SDFAI calculation and plot:

```bash
python analysis/00_runtime/01_run_analysis.py 04.01
python analysis/00_runtime/01_run_analysis.py 04.02
python analysis/00_runtime/01_run_analysis.py 04.03
python analysis/00_runtime/01_run_analysis.py 04.04 -- --boxed
python analysis/00_runtime/01_run_analysis.py 04.06 -- /path/to/monthly_environment.xlsx /path/to/environment_outputs

CCF_SUPPLEMENT_MODE=TOC_WT python analysis/00_runtime/01_run_analysis.py 20.02
CCF_SUPPLEMENT_MODE=AIR_TEMP python analysis/00_runtime/01_run_analysis.py 20.02
CCF_SUPPLEMENT_MODE=N_P_RATIO python analysis/00_runtime/01_run_analysis.py 20.02

python analysis/00_runtime/01_run_analysis.py 22.01 -- WL formal
python analysis/00_runtime/01_run_analysis.py 22.04 -- --unambiguous
python analysis/00_runtime/01_run_analysis.py 15.02 -- /path/to/overall_inputs /path/to/candidate_inputs /path/to/regression_outputs
```

CCF plots consume `level_sensitivity` results. Formal TITAN2 uses 500 permutations and 500 bootstraps. Run the threshold and ridgeline plots before the label and alignment exports. Plotting scripts require the tables and cached predictions produced by upstream calculations.

For `15.02`, both input folders require `groupsize.csv` and `env.csv`. The overall folder also requires `alpha1.csv`; the candidate folder requires `04_pathogen_candidate_abundance.csv`. The complete `tpm_rarefied.csv` must be two directory levels above the candidate folder. This candidate-richness calculation uses month-clustered HC1 standard errors. Use the complete-community outputs of `15.01` and the candidate outputs of `15.02`.

`21.02` accepts `input_directory output_directory`. `21.04` accepts `alpha_csv group_csv fraction_abundance_csv figure_directory table_directory`. The composition function set uses the rarefied candidate subset by default; `23.01` exports candidate genus composition from original counts.
