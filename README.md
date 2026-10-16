# Danjiangkou 16S analysis

The `analysis` directory contains the manuscript analyses in execution order. Each script includes its associated calculations and plots where they share inputs. `gamm_functions.R` provides shared model functions; `run.py` is the single launcher.

| Step | Directory | Analysis |
| --- | --- | --- |
| 01 | [sequence_processing](analysis/01_sequence_processing/) | ASV inference, taxonomy, rooted tree and MBPD candidate screening |
| 02 | [data_preparation](analysis/02_data_preparation/) | SILVA tables, whole-community rarefaction and candidate extraction |
| 03 | [hydrology_environment](analysis/03_hydrology_environment/) | Discharge, Gamma-MLE SRI/SDFAI, water level, environment and sampling map |
| 04 | [community_composition](analysis/04_community_composition/) | Shared ASVs, taxonomy network and taxonomic composition |
| 05 | [alpha_diversity](analysis/05_alpha_diversity/) | Diversity comparisons and richness trends |
| 06 | [beta_diversity](analysis/06_beta_diversity/) | PCoA comparisons and temporal trends |
| 07 | [environmental_associations](analysis/07_environmental_associations/) | LEfSe, distance decay, RDA, Mantel, variation partitioning and random forest |
| 08 | [water_level_regression](analysis/08_water_level_regression/) | Community richness, candidate abundance and candidate richness regressions |
| 09 | [community_assembly](analysis/09_community_assembly/) | iCAMP and assembly-process patterns |
| 10 | [networks](analysis/10_networks/) | Static networks, node removal and temporal networks |
| 11 | [cyanobacteria_ccf](analysis/11_cyanobacteria_ccf/) | Cyanobacterial dynamics and cross-correlations |
| 12 | [candidate_dynamics](analysis/12_candidate_dynamics/) | Candidate abundance, richness and community fraction |
| 13 | [indicator_thresholds](analysis/13_indicator_thresholds/) | TITAN2, threshold curves and ridgelines |

## Requirements and inputs

Use R 4.5.1 and Python 3.9 or later. R package versions are listed in [R_packages.csv](analysis/R_packages.csv); Python dependencies are in [requirements.txt](analysis/requirements.txt). Sequence processing requires Linux, Bash and QIIME 2; candidate screening requires NCBI BLAST+.

Keep the study dataset in a separate directory with its existing relative paths and column names. ASV tables have ASVs in rows and samples in columns. Sample identifiers must match the group and environment tables. Original counts are used for shared-ASV, taxonomy-network and iCAMP analyses. Diversity analyses use the complete community rarefied once to 21,000 reads per sample; candidate richness uses a subset of that same table. Environment inputs must already use the units specified by their column labels. Reference databases and input data are supplied separately.

## Run

Run commands from the repository root. Analysis IDs and script paths are listed in [workflow.tsv](analysis/workflow.tsv).

```bash
python analysis/run.py --list
python analysis/run.py --list --step 03
python analysis/run.py 02.01 --data-root /path/to/study_data --check
python analysis/run.py 02.01 --data-root /path/to/study_data
```

`--check` checks the script and data-directory paths. Use `--rscript /path/to/Rscript` or `--bash /path/to/bash` to select an executable. Arguments after `--` are passed to the analysis. Outputs are written within the study data directory. The launcher sets `DJK_CODE_ROOT` and `DJK_PROJECT_ROOT`.

For repeated runs, set `DJK_PROJECT_ROOT` once (`export DJK_PROJECT_ROOT=/path/to/study_data` in Bash, or `$env:DJK_PROJECT_ROOT = 'C:/path/to/study_data'` in PowerShell).

## Sequence processing and assembly

The sequence workflow expects paired `sample_R1.fastq.gz` and `sample_R2.fastq.gz` files. Set `ENV_NAME` and `CONDA_SH` for the QIIME installation and use a new output directory.

```bash
CLEAN_DIR=/path/to/primer_trimmed_reads \
BASE_WORK=/path/to/results/qiime2 \
RUN_TAG=analysis \
CLASSIFIER_QZA=/path/to/silva_classifier.qza \
bash analysis/01_sequence_processing/01_process_amplicons.sh

python analysis/01_sequence_processing/02_screen_candidates.py \
  --asv-fasta /path/to/rep_seqs.fasta \
  --db-fasta /path/to/MBPD_reference.fasta \
  --db-tax /path/to/MBPD_taxonomy.tsv \
  --outdir /path/to/results/mbpd \
  --identity 97 --coverage 95 --threads 8 --max-target-seqs 50

ASV_EXPORT_DIR=/path/to/06_export \
ICAMP_OUTPUT_DIR=/path/to/results/icamp \
Rscript analysis/09_community_assembly/01_icamp.R
```

MBPD requires `--coverage 95` for the manuscript filter; its default is 90. The E-value threshold is `1e-20`. `--abundance` optionally supplies an ASV count CSV. iCAMP uses seed 123, 1,000 randomizations, minimum bin size 48, 16 workers and a 50-GB memory allowance. `ASV_TABLE_FILE` and `TREE_FILE` override its individual inputs.

## Analysis-specific arguments

The following commands assume `DJK_PROJECT_ROOT` is set:

```bash
python analysis/run.py 03.01
python analysis/run.py 03.02
python analysis/run.py 03.04 -- /path/to/monthly_environment.xlsx /path/to/environment_outputs
python analysis/run.py 08.01 -- /path/to/overall_inputs /path/to/candidate_inputs /path/to/regression_outputs
python analysis/run.py 08.02 -- /path/to/candidate_inputs /path/to/candidate_richness_outputs
python analysis/run.py 11.03
python analysis/run.py 13.01 -- WL formal
python analysis/run.py 13.02
```

`03.02` calculates calendar-month Gamma-MLE SRI and adjacent-month SDFAI and creates the final plot. `07.06` fits the final 19-predictor random forest excluding Q and plots its results. `08.01` requires `alpha1.csv`, `groupsize.csv` and `env.csv` in its overall folder, plus candidate counts and matching metadata in its candidate folder. `08.02` uses month-clustered HC1 standard errors; its candidate folder contains `groupsize.csv`, `env.csv` and `04_pathogen_candidate_abundance.csv`, with the complete-community `tpm_rarefied.csv` two directory levels above it.

`11.03` calculates TOC/WT, air-temperature and TN:TP cross-correlations and their plots. `12.01` accepts `input_directory output_directory`; `12.02` accepts `alpha_csv group_csv fraction_abundance_csv figure_directory table_directory`. `13.01` in formal mode uses 500 permutations and 500 bootstraps; `13.02` creates the threshold and ridgeline plots with the final genus labels. Plotting requires upstream results and cached predictions. `04.04` also exports candidate genus composition from original counts.
