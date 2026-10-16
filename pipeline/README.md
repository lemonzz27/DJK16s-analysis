# Sequence processing and community assembly

These scripts accompany *Size dependent microbial responses and risk signals during a reservoir drought to flood transition*. They prepare ASV inputs, screen candidate pathogen-related sequences and estimate community assembly. Downstream statistical analysis and plotting scripts are provided elsewhere in this repository.

## Requirements

Run the shell workflows on Linux with Bash 4.4 or later. They target the QIIME 2/q2cli 2026.4 amplicon distribution and require the corresponding DADA2, feature-classifier, alignment and phylogeny plugins, BIOM tools, Python, scikit-learn, gzip and GNU command-line utilities. Use the intended QIIME 2 environment; the scripts do not install or update it. Set `ENV_NAME` if its name differs from `rachis-qiime2-2026.4`, or `CONDA_SH` if Conda initialization is outside its usual location.

MBPD screening requires Python 3.9 or later and NCBI BLAST+ (`makeblastdb` and `blastn`). It uses only the Python standard library. The iCAMP script requires R with `iCAMP` and `ape`; its retained configuration uses 16 workers and a 50-GB memory allowance.

The installed versions and complete reference database files used in the original analysis have not been independently reconstructed. QIIME 2 records its environment when run, and the iCAMP script records `sessionInfo()`. Retain these records with the outputs.

## Input preparation

Provide primer-trimmed paired-end reads named `sample_R1.fastq.gz` and `sample_R2.fastq.gz`, either directly under the input directory or one directory level below it. `MAX_INPUT_DEPTH` can be increased for a deeper directory layout. Sample IDs must be unique. Optional `METADATA_FILE` must be a QIIME 2 metadata TSV whose first column contains matching sample IDs.

Taxonomy uses SILVA 138.2 SSU NR99 references for the V3-V4 region. Provide either a compatible trained classifier or the matching reference sequence and taxonomy artifacts needed to train it. Reference data are external inputs; download provenance, release identifiers and file checksums must be retained. These scripts do not create reference data from an unspecified database release.

## Read processing

Run from the repository root, replacing example input paths with your own:

```bash
CLEAN_DIR=/path/to/primer_trimmed_reads \
BASE_WORK="$PWD/results/qiime2" \
RUN_TAG=manuscript \
CLASSIFIER_QZA=/path/to/silva_138.2_v3v4_classifier.qza \
AMPLICON_REGION=V3-V4 \
TAXONOMY_DB='SILVA 138.2 SSU NR99' \
bash pipeline/16S.sh
```

The output is `results/qiime2/ASV_QIIME2_manuscript`. The script stops if its `WORK` directory already exists. Choose a new `RUN_TAG` or `WORK` for another run.

Default DADA2 settings are no left trimming and no fixed-length truncation, maximum expected errors of 2 per read, quality truncation threshold 2, minimum overlap 12 nt, zero merge mismatches, independent pooling, 1,000,000 reads for error learning and consensus chimera removal. These settings assume the input reads have already had primers removed. Thread counts and settings are environment variables; changes should be recorded with the resulting analysis.

Classification defaults are confidence 0.7 and automatic read orientation. Tree inference uses MAFFT, masking with maximum gap frequency 1.0 and minimum conservation 0.4, FastTree with one thread, and midpoint rooting. `RUN_TAXONOMY=0` and `RUN_TREE=0` can be used to stop after ASV processing and exports when the next workflow will supply taxonomy and phylogeny.

Main outputs include `06_export/ASV_table.tsv`, `06_export/rep_seqs/dna-sequences.fasta`, taxonomy tables and `06_export/phylogeny/rooted_tree/tree.nwk` when enabled. `00_provenance` contains parameters, environment information, checksums and the executed script.

## Continue from existing ASV artifacts

This alternative starts from `02_dada2` and `03_qc` of an existing QIIME 2 run. It can train a SILVA classifier using the current scikit-learn version.

```bash
WORK="$PWD/results/qiime2/ASV_QIIME2_manuscript" \
REF_SEQS=/path/to/silva-138.2-ssu-nr99-seqs-V3-V4-uniq.qza \
REF_TAX=/path/to/silva-138.2-ssu-nr99-tax-V3-V4-derep-uniq.qza \
RESULT_DIR="$PWD/results/taxonomy_tree_manuscript" \
bash pipeline/continue_ASV_taxonomy_tree.sh
```

`WORK` must contain the ASV table, representative sequences, denoising statistics and base-transition statistics under `02_dada2`, and ASV/sample frequency artifacts under `03_qc`. Outputs are written to the new `RESULT_DIR`, which defaults to `WORK/taxonomy_tree`. Its `06_export` directory is a complete export suitable for the following steps. Existing `RESULT_DIR` paths are refused. The input run is not modified or cleaned.

`SILVA_DIR` can be used instead of explicit reference artifact paths when its filenames match the example. `CLASSIFIER_DIR` can point to a classifier cache; it defaults to `RESULT_DIR/reference_classifier`. A cached classifier that fails validation causes the script to stop. No files are deleted automatically. Classification and tree defaults match the read-processing workflow.

## Candidate pathogen-related sequences

```bash
python pipeline/mbpd_asv_direct.py \
  --asv-fasta results/taxonomy_tree_manuscript/06_export/rep_seqs/dna-sequences.fasta \
  --db-fasta /path/to/MBPD_reference.fasta \
  --db-tax /path/to/MBPD_reference_taxonomy.tsv \
  --outdir results/mbpd_manuscript \
  --identity 97 \
  --coverage 95 \
  --threads 8 \
  --max-target-seqs 50
```

The example uses the manuscript's 97% identity, 95% coverage and fixed E-value threshold of `1e-20`. The original script's default coverage is retained at 90%, so `--coverage 95` must be passed explicitly for this example. The available filtered-hit table contains 57,964 hits for 2,755 ASVs: its minimum recorded qcovs is 95%, minimum identity is 97%, and maximum E-value is 1.2e-150. This supports the recorded downstream set. The historical BLAST command and HSP-level coverage cannot be reconstructed from that filtered table alone.

The search uses megablast on both strands, with a maximum of 50 target sequences. Output includes all returned alignments, annotated alignments and a best-score summary that retains ambiguous annotations. Coverage is applied through BLAST's `-qcov_hsp_perc` option; the output also records `qcovs`. This script does not implement lowest-common-ancestor taxonomy assignment. Any downstream LCA taxonomy table needs its own documented preparation step. Sequence matches are candidate pathogen-related signals and do not establish viability or pathogenicity.

The output directory must be new or empty. `--abundance` optionally accepts a comma-delimited count table with ASV IDs in its first column; it is not a TSV reader. The optional abundance export retains the input header and candidate rows without rarefaction or normalization.

## Community assembly

```bash
ASV_EXPORT_DIR="$PWD/results/taxonomy_tree_manuscript/06_export" \
ICAMP_OUTPUT_DIR="$PWD/results/icamp_manuscript" \
Rscript pipeline/run_icamp_bin48.R
```

`ASV_TABLE_FILE` and `TREE_FILE` can override the two input paths individually. `ICAMP_R_LIBRARY` can specify an additional R package library. `ICAMP_OUTPUT_DIR` defaults to `results/icamp_bin48_rand1000` relative to the current directory and must be new or empty.

The workflow reads the original ASV count table, removes zero-abundance samples/taxa, matches taxa to tree tips and saves the exact community and tree used. It does not rarefy the data. Scientific settings are retained from the source: seed 123, 1,000 randomizations, minimum bin size 48, `ds = 0.2`, `phylo.metric = "bMPD"`, phylogenetic randomization within bins, taxonomic randomization across all taxa, `sig.index = "SES.RC"`, SES cutoff 1.96 and RC cutoff 0.95. All other `icamp.big` arguments are preserved. Outputs include the complete R result, pairwise process importance, the overall mean process summary, input QC, parameters and the R session information.

## Validation scope

The publication copies were checked with Bash syntax validation, R parsing, Python compilation and focused checks of path guards and the FASTA reader. The iCAMP call and BLAST search function were compared with their source versions. Full QIIME 2, BLAST and iCAMP computations were not rerun during code preparation. These checks establish source validity and the stated limited behavior; they do not establish end-to-end reproduction of the manuscript results.
