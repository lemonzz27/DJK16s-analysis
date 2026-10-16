# Run guide

Use a separate working copy with the study inputs staged at the paths in `provenance/input_files.tsv`. Several original scripts write into their input workspaces and can replace generated results. The input catalogue is not an exhaustive list of dynamically named caches.

## Execution

Set an absolute `DJK_PROJECT_ROOT` and ensure Rscript is on PATH (or pass `--rscript`). The launcher sets both project-root variables, selects the recorded working directory and starts a separate R process. Its use of `source()` matches the original execution method for scripts that infer paths from the working directory. Do not directly call those scripts with `Rscript path/to/script.R`, because `--file` can make their source directory become their input directory.

```bash
python scripts/run_analysis.py --list
python scripts/run_analysis.py A001 --check
python scripts/run_analysis.py A001
```

`--check` reports the selected script, working directory and whether they exist; it does not validate every input or run the analysis. An unavailable data workspace causes a clear error. Use `--working-directory` only when the corresponding input set has been deliberately staged there. For scripts that require positional options, place their options after `--`.

For Python source-preparation/layout scripts, run the script directly from the repository root after staging its inputs. For QIIME 2, MBPD and iCAMP use `pipeline/README.md`.

## Manuscript execution order

1. Sequence processing, taxonomy/tree construction and MBPD screening precede the single whole-community rarefaction. Prepare sample groups, taxonomy and the updated environment table before ecological analyses. The two taxonomy sources supported by legacy helpers are not interchangeable with the manuscript's SILVA analysis.
2. Figure 1A uses `prepare_and_verify.py prepare`, then `calculate_indices.R` to form the complete monthly grid, `baihe_sri_matrix_20260922/calculate.R`, and `baihe_gamma_sdfai_20260922/export_plot.R --boxed`. The selected manuscript result is Gamma MLE / SDFAI; the earlier preliminary SRI columns are not final results. The final analysis has 510 observed months, 506 adjacent pairs, 19 dry-to-wet and 21 wet-to-dry transitions. Figure 1B and C-F use the water-level and monthly-environment plotting scripts, including the corrected Chla units.
3. Figure 2 uses original counts for shared/unique ASVs and the radial network; monthly composition and richness use the single 21,000-read rarefaction. Run alpha/PCoA calculations before their time-trend and figure-export scripts.
4. Figure 3 uses `work/rda_cca_noQ/run_rda_cca_noQ.R` with its staged 19-variable input, followed by the full-arrow RDA plotting script. Run distance-decay and Mantel with their corresponding 19-variable working inputs. The older 24-variable environment-update entry point is not the final manuscript analysis.
5. Figure 4 starts with `pipeline/run_icamp_bin48.R`, then process aggregation and GAMM plotting. Static and moving-window network calculations precede the large-font exports. Their outputs/caches must be staged at the paths expected by the figure helpers.
6. Figure 5A uses only cyanobacterial extraction and monthly plotting from the directory named SEM; no SEM fit is part of the current manuscript. CCF uses the `28_...CCF` scripts. The supplement calculation must be run with `CCF_SUPPLEMENT_MODE=TOC_WT`, `AIR_TEMP` and `N_P_RATIO` to generate the separate relations. The main manuscript uses the original-level sensitivity outputs, not the alternative first-difference outputs. The final Python layout helper consumes already exported PDF panels.
7. Figure 6 uses candidate GAMMs and the formal ASV-level TITAN2 run, followed by sum(z), ridgeline, alignment and label exports. TITAN2 retains 476 ASVs and reports 181 robust indicators from 2,755 candidates; the formal settings use 500 permutations and 500 bootstraps. Run the final label script with `--unambiguous`. Its code-fragment dependencies were checked after reformatting.
8. Supplementary Fig. 10 uses the October 2 random-forest rerun and plotter with 19 predictors excluding Q. Supplementary Fig. 15 uses `S15_candidate_richness_single_rarefaction.R` with three positional paths: overall regression inputs, candidate regression inputs and output directory. Candidate richness comes from the complete-community rarefaction and uses month-clustered HC1 standard errors. The VPA label correction is also included.

Final PowerPoint/manual multi-panel assembly is distinct from numerical and single-panel reproduction. A single automatic command that regenerates the complete submitted manuscript has not been established.

## Selected scripts and working directories

| ID | Analysis | Working directory |
| --- | --- | --- |
| A001 | Single rarefaction and candidate subset | `.` |
| A002 | Taxonomy conversion | `.` |
| A003 | Alpha diversity calculation | `99_运行环境/串行运行工作区/SILVA/a多样性` |
| A004 | Supplementary Fig 5 alpha comparison | `99_运行环境/串行运行工作区/SILVA/a多样性` |
| A005 | Shared GAMM functions | `.` |
| A006 | Fig 2C richness GAMM calculation | `99_运行环境/串行运行工作区/SILVA/a多样性年际变化` |
| A007 | Water-level GAMM inset | `99_运行环境/串行运行工作区/SILVA/a多样性年际变化` |
| A008 | Fig 2C richness presentation | `.` |
| A009 | Fig 6A candidate time trends | `.` |
| A010 | Fig 2D month and size PCoA | `.` |
| A011 | Supplementary Fig 6 site PCoA | `99_运行环境/串行运行工作区/SILVA/PCOA` |
| A012 | Supplementary Fig 13 candidate PCoA | `99_运行环境/串行运行工作区/SILVA/致病菌/PCOA` |
| A013 | PCoA1 time course and RF response preparation | `99_运行环境/串行运行工作区/SILVA/PCOA1年际变化` |
| A014 | PCoA GAMM helper bridge | `.` |
| A015 | Fig 2A complete-community Venn | `.` |
| A016 | Supplementary Fig 12 candidate Venn | `.` |
| A017 | Fig 2A radial taxonomy network | `.` |
| A018 | Supplementary Fig 11 candidate relative abundance | `.` |
| A019 | MBPD grouped composition functions | `.` |
| A020 | Fig 2B and Supplementary Figs 3-4 monthly composition | `.` |
| A021 | Fig 4F-G static networks and node removal | `99_运行环境/串行运行工作区/SILVA/network` |
| A022 | Fig 4G node-removal figure | `.` |
| A023 | Static network topology figure | `.` |
| A024 | Static-network legend | `.` |
| A025 | Fig 4C-E sliding-window networks | `99_运行环境/串行运行工作区/SILVA/network年际变化` |
| A026 | Network GAMM helper bridge | `.` |
| A027 | Fig 4C-E final large-axis figures | `.` |
| A028 | Supplementary Fig 2 LEfSe | `99_运行环境/串行运行工作区/SILVA/lefse` |
| A029 | Fig 3D-F RDA calculation | `.` |
| A030 | Fig 3D-F final full-arrow RDA plotting | `.` |
| A031 | Fig 3G Mantel and environment correlations | `.` |
| A032 | Fig 3A-C environmental/spatial/temporal decay | `99_运行环境/串行运行工作区/SILVA/decay` |
| A033 | Fig 3A-C combined plot | `.` |
| A034 | Supplementary Fig 8 water-level decay | `.` |
| A035 | Supplementary Fig 10 / Data table 15 RF | `.` |
| A036 | Supplementary Fig 10 RF plotting | `.` |
| A037 | RF predecessor provenance | `99_运行环境/串行运行工作区/SILVA/随机森林` |
| A038 | Supplementary Figs 7 and 16 sample-level water regression | `.` |
| A039 | Supplementary Fig 15 / Data table 26 candidate richness regression | `.` |
| A040 | Supplementary Fig 9 / Data table 14 VPA | `.` |
| A041 | Supplementary Fig 9 corrected VPA labels | `.` |
| A042 | Fig 4A iCAMP aggregation | `.` |
| A043 | Fig 4B stochastic/deterministic GAMM | `.` |
| A044 | Fig 4A-B iCAMP large-font exports | `.` |
| A045 | Fig 4A final frame and legend | `.` |
| A046 | Supplementary Data full composition and alpha exports | `.` |
| A047 | Original radial-network calculation helper | `.` |
| A048 | Fig 6A cached-layout source | `.` |
| A049 | Fig 6A wider plotting helper | `.` |
| A050 | Fig 6A-B final aligned panels | `.` |
| A051 | Figure 1A: Baihe discharge source preparation | `.` |
| A052 | Figure 1A: complete monthly grid preparation | `.` |
| A053 | Figure 1A: calendar-month Gamma MLE SRI-1 and SDFAI | `.` |
| A054 | Figure 1A: final SDFAI export | `.` |
| A055 | Figure 1C-F: monthly HRT, TP, TN and chlorophyll a | `.` |
| A056 | Figure 1B: monthly water level and net inflow | `.` |
| A057 | Figure 1F source conversion | `.` |
| A058 | Figure 1F vector export | `.` |
| A059 | Figure 5B,C,E: water level, Chla and cyanobacteria CCF | `.` |
| A060 | Figure 5D,F,G: TN:TP and temperature CCF | `.` |
| A061 | Figure 5 original-level CCF panel plotting | `.` |
| A062 | Figure 5 TN:TP and temperature panel plotting | `.` |
| A063 | Cyanobacteria extraction configuration | `.` |
| A064 | Figure 5A: extract sample-level cyanobacteria | `.` |
| A065 | Figure 5A: monthly cyanobacterial abundance | `.` |
| A066 | Figure 5A local launcher | `.` |
| A067 | Figure 5 final vector-panel layout | `.` |
| A068 | Figure 6B,C: ASV-level TITAN2 | `.` |
| A069 | Figure 6B: ASV community sum(z) | `.` |
| A070 | Figure 6C: top 20 unique-species representative ASVs | `.` |
| A071 | Figure 6C: final unambiguous genus labels | `.` |
| A072 | ASV TITAN2 numerical QA | `.` |
| A073 | Local TITAN ridge launcher | `.` |
| A074 | Supplementary Fig. 11: full richness and candidate fraction | `.` |
| A075 | Supplementary Fig. 1: four-station study map | `.` |
| A076 | GIS source inspection | `.` |
