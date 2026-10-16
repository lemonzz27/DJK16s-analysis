# Danjiangkou Reservoir 16S analysis

Analysis code for **Size dependent microbial responses and risk signals during a reservoir drought to flood transition**.

Repository: [https://github.com/lemonzz27/DJK16s-analysis](https://github.com/lemonzz27/DJK16s-analysis) (private).

Version: `0.1.0-review`, prepared on 16 October 2026 against the Nature Communications submission manuscript and Supplementary Information of the same date. This is a private code preparation version. It is not a public archival release and has no DOI.

The repository contains 80 selected analysis, plotting and sequence-processing scripts: 74 R scripts, four Python scripts and two Bash scripts. Selection used manuscript methods, reported values, saved results and source dependencies, rather than modification dates alone. Original source files were left unchanged. The selection evidence and original SHA-256 checksums are in [the provenance manifest](provenance/source_selection.tsv).

## Contents

- [Analysis index](docs/analysis_index.md): script-to-analysis and figure mapping.
- [Run guide](docs/run_guide.md): execution order, working directories and important settings.
- [Sequence processing](pipeline/README.md): QIIME 2, MBPD and iCAMP inputs and commands.
- [Input requirements](docs/inputs.md): required data, cached results and directory layout.
- [Reproduction notes](docs/reproduction_notes.md): script-specific limitations and provenance details.
- [Validation](docs/validation.md): checks completed during preparation and their limits.
- [Release requirements](docs/release_requirements.md): licence, reviewer access and archival version.

## Environment

The local R environment checked during preparation is R 4.5.1. [R_packages.csv](environment/R_packages.csv) records observed installed versions, not a reconstructed historical lockfile. Python support scripts use openpyxl, NumPy and PyMuPDF; MBPD screening uses the Python standard library plus NCBI BLAST+. The shell workflows target Linux and QIIME 2 2026.4. See [environment notes](environment/README.md).

Set `DJK_PROJECT_ROOT` to the absolute path of this checkout. Paths that previously pointed to the author's project have been replaced by this setting. If it is unset, scripts use `.`. Some older analyses read filenames relative to their working directory: use the working directories in the run guide rather than running every script from the same folder.

The source project layout and names are retained because scripts share data and helper files through these paths. `work`, `work_artifact` and `outputs` contain selected scripts required by the final manuscript, including later corrections; these names do not designate discarded versions in this repository.

## Scientific scope

The selected versions include the Gamma-MLE Baihe SRI/SDFAI calculation, the October 2026 19-predictor random forest rerun without Q, the corrected single-rarefaction candidate-richness regression, ASV-level TITAN2 and final figure export helpers. Exploratory DLNM, MIP, BugBase and superseded genus/species TITAN analyses are not manuscript entry points and are excluded.

Source comments and decorative section banners were removed; interpreter directives were retained. R code was formatted consistently, and parsing was used to verify that calculations were preserved after the documented path substitutions. Pipeline changes additionally parameterize inputs and prevent reuse or automatic deletion of existing output directories. Details are in [pipeline changes](docs/pipeline_changes.md).

The filtered MBPD hit table was checked: 57,964 hits cover 2,755 ASVs; all recorded query coverages are at least 95% and identities at least 97%. The script's original 90% default was retained; the manuscript command explicitly uses `--coverage 95`. The historical BLAST command and HSP-level coverage still require the original run record.

Input data, reference databases and analysis caches are not included in this code-only repository. Full sequence processing and all manuscript analyses were not rerun during packaging. Some final multi-panel layouts include manual assembly; see the run guide for the reproducible statistical and plotting components.

## Licence and citation

No new open-source licence has been granted for this private preparation version. The authors must select a licence and confirm software authorship before public release. Third-party packages and reference datasets retain their own terms. See [citation information](CITATION.md).
