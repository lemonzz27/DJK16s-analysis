#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

ENV_NAME="${ENV_NAME:-rachis-qiime2-2026.4}"
EXPECTED_QIIME_VERSION="${EXPECTED_QIIME_VERSION:-2026.4}"
CONDA_SH="${CONDA_SH:-}"

CLEAN_DIR="${CLEAN_DIR:?Set CLEAN_DIR to the primer-trimmed paired-end FASTQ directory}"
BASE_WORK="${BASE_WORK:-${DJK_PROJECT_ROOT:-${PWD}}/results/qiime2}"
METADATA_FILE="${METADATA_FILE:-}"

CLASSIFIER_QZA="${CLASSIFIER_QZA:-}"
AMPLICON_REGION="${AMPLICON_REGION:-unknown}"
TAXONOMY_DB="${TAXONOMY_DB:-unknown}"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
TAX_N_JOBS="${TAX_N_JOBS:-4}"
MAFFT_THREADS="${MAFFT_THREADS:-${THREADS}}"
FASTTREE_THREADS="${FASTTREE_THREADS:-1}"

TRIM_LEFT_F="${TRIM_LEFT_F:-0}"
TRIM_LEFT_R="${TRIM_LEFT_R:-0}"

TRUNC_F="${TRUNC_F:-0}"
TRUNC_R="${TRUNC_R:-0}"
MAX_EE_F="${MAX_EE_F:-2}"
MAX_EE_R="${MAX_EE_R:-2}"
TRUNC_Q="${TRUNC_Q:-2}"
MIN_OVERLAP="${MIN_OVERLAP:-12}"
MAX_MERGE_MISMATCH="${MAX_MERGE_MISMATCH:-0}"
POOLING_METHOD="${POOLING_METHOD:-independent}"
N_READS_LEARN="${N_READS_LEARN:-1000000}"

TAX_CONFIDENCE="${TAX_CONFIDENCE:-0.7}"
TAX_READ_ORIENTATION="${TAX_READ_ORIENTATION:-auto}"
MASK_MAX_GAP_FREQUENCY="${MASK_MAX_GAP_FREQUENCY:-1.0}"
MASK_MIN_CONSERVATION="${MASK_MIN_CONSERVATION:-0.4}"

RUN_TAXONOMY="${RUN_TAXONOMY:-1}"
RUN_TREE="${RUN_TREE:-1}"
CHECK_GZIP="${CHECK_GZIP:-1}"
MAX_INPUT_DEPTH="${MAX_INPUT_DEPTH:-2}"

RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)_$$}"
WORK="${WORK:-${BASE_WORK}/ASV_QIIME2_${RUN_TAG}}"

die() {
    printf '[错误] %s\n' "$*" >&2
    exit 1
}

warn() {
    printf '[警告] %s\n' "$*" >&2
}

require_positive_integer() {
    local name="$1"
    local value="$2"
    [[ "${value}" =~ ^[1-9][0-9]*$ ]] || die "${name} 必须为正整数，当前值：${value}"
}

require_nonnegative_integer() {
    local name="$1"
    local value="$2"
    [[ "${value}" =~ ^[0-9]+$ ]] || die "${name} 必须为非负整数，当前值：${value}"
}

require_nonnegative_number() {
    local name="$1"
    local value="$2"
    [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] || \
        die "${name} 必须为非负数，当前值：${value}"
}

require_unit_interval() {
    local name="$1"
    local value="$2"
    require_nonnegative_number "${name}" "${value}"
    awk -v value="${value}" 'BEGIN {exit !(value >= 0 && value <= 1)}' || \
        die "${name} 必须在 0 到 1 之间，当前值：${value}"
}

make_absolute() {
    local path="$1"
    if [[ "${path}" == /* ]]; then
        printf '%s\n' "${path}"
    else
        printf '%s/%s\n' "${START_DIR}" "${path}"
    fi
}

activate_qiime_environment() {
    if [[ "${CONDA_DEFAULT_ENV:-}" == "${ENV_NAME}" ]]; then
        return 0
    fi

    local conda_init="${CONDA_SH}"
    if [[ -z "${conda_init}" ]] && command -v conda >/dev/null 2>&1; then
        conda_init="$(conda info --base)/etc/profile.d/conda.sh"
    fi

    if [[ -z "${conda_init}" ]]; then
        local candidate
        for candidate in \
            "${HOME}/miniconda3/etc/profile.d/conda.sh" \
            "${HOME}/anaconda3/etc/profile.d/conda.sh" \
            "/opt/conda/etc/profile.d/conda.sh"; do
            if [[ -f "${candidate}" ]]; then
                conda_init="${candidate}"
                break
            fi
        done
    fi

    [[ -n "${conda_init}" && -f "${conda_init}" ]] || \
        die "找不到 conda.sh；请用 CONDA_SH=/绝对路径/conda.sh 指定。"

    source "${conda_init}"
    conda activate "${ENV_NAME}"
}

on_error() {
    local exit_code="$1"
    local line_no="$2"
    local command_text="$3"
    printf '[失败] exit=%s line=%s command=%s\n' \
        "${exit_code}" "${line_no}" "${command_text}" >&2
}

on_exit() {
    local exit_code=$?
    if [[ "${exit_code}" -eq 0 ]]; then
        printf '[完成] %s\n' "$(date '+%F %T')"
    else
        printf '[终止] exit=%s；日志：%s\n' "${exit_code}" "${LOG_FILE:-尚未创建}" >&2
    fi
}

START_DIR="$(pwd -P)"
SCRIPT_SOURCE="${BASH_SOURCE[0]}"
if [[ "${SCRIPT_SOURCE}" != /* ]]; then
    if [[ -f "${START_DIR}/${SCRIPT_SOURCE}" ]]; then
        SCRIPT_SOURCE="${START_DIR}/${SCRIPT_SOURCE}"
    else
        RESOLVED_SCRIPT="$(command -v "${SCRIPT_SOURCE}" || true)"
        [[ -n "${RESOLVED_SCRIPT}" ]] || die "无法定位当前脚本文件：${SCRIPT_SOURCE}"
        if [[ "${RESOLVED_SCRIPT}" == /* ]]; then
            SCRIPT_SOURCE="${RESOLVED_SCRIPT}"
        else
            SCRIPT_SOURCE="${START_DIR}/${RESOLVED_SCRIPT}"
        fi
    fi
fi
CLEAN_DIR="$(make_absolute "${CLEAN_DIR}")"
BASE_WORK="$(make_absolute "${BASE_WORK}")"
WORK="$(make_absolute "${WORK}")"
if [[ -n "${METADATA_FILE}" ]]; then
    METADATA_FILE="$(make_absolute "${METADATA_FILE}")"
fi
if [[ -n "${CLASSIFIER_QZA}" ]]; then
    CLASSIFIER_QZA="$(make_absolute "${CLASSIFIER_QZA}")"
fi

for flag_name in RUN_TAXONOMY RUN_TREE CHECK_GZIP; do
    flag_value="${!flag_name}"
    [[ "${flag_value}" == 0 || "${flag_value}" == 1 ]] || \
        die "${flag_name} 只能设为 0 或 1。"
done

require_positive_integer THREADS "${THREADS}"
require_positive_integer TAX_N_JOBS "${TAX_N_JOBS}"
require_positive_integer MAFFT_THREADS "${MAFFT_THREADS}"
require_positive_integer FASTTREE_THREADS "${FASTTREE_THREADS}"
require_positive_integer MAX_INPUT_DEPTH "${MAX_INPUT_DEPTH}"
require_nonnegative_integer TRIM_LEFT_F "${TRIM_LEFT_F}"
require_nonnegative_integer TRIM_LEFT_R "${TRIM_LEFT_R}"
require_nonnegative_integer TRUNC_F "${TRUNC_F}"
require_nonnegative_integer TRUNC_R "${TRUNC_R}"
require_nonnegative_integer TRUNC_Q "${TRUNC_Q}"
require_nonnegative_number MAX_EE_F "${MAX_EE_F}"
require_nonnegative_number MAX_EE_R "${MAX_EE_R}"
require_positive_integer MIN_OVERLAP "${MIN_OVERLAP}"
require_nonnegative_integer MAX_MERGE_MISMATCH "${MAX_MERGE_MISMATCH}"
require_positive_integer N_READS_LEARN "${N_READS_LEARN}"
if [[ "${TAX_CONFIDENCE}" != disable ]]; then
    require_unit_interval TAX_CONFIDENCE "${TAX_CONFIDENCE}"
fi
require_unit_interval MASK_MAX_GAP_FREQUENCY "${MASK_MAX_GAP_FREQUENCY}"
require_unit_interval MASK_MIN_CONSERVATION "${MASK_MIN_CONSERVATION}"

[[ "${MIN_OVERLAP}" -ge 4 ]] || die "QIIME 2 要求 MIN_OVERLAP >= 4。"
[[ "${POOLING_METHOD}" == independent || "${POOLING_METHOD}" == pseudo ]] || \
    die "POOLING_METHOD 只能是 independent 或 pseudo。"
case "${TAX_READ_ORIENTATION}" in
    same|reverse-complement|auto|both) ;;
    *) die "TAX_READ_ORIENTATION 必须是 same、reverse-complement、auto 或 both。" ;;
esac

[[ -d "${CLEAN_DIR}" && -r "${CLEAN_DIR}" ]] || \
    die "CleanData 目录不存在或不可读：${CLEAN_DIR}"
if [[ -n "${METADATA_FILE}" ]]; then
    [[ -f "${METADATA_FILE}" && -r "${METADATA_FILE}" ]] || \
        die "metadata 文件不存在或不可读：${METADATA_FILE}"
fi
if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
    [[ -n "${CLASSIFIER_QZA}" ]] || \
        die "RUN_TAXONOMY=1 时必须设置 CLASSIFIER_QZA=/绝对路径/classifier.qza。"
    [[ -f "${CLASSIFIER_QZA}" && -r "${CLASSIFIER_QZA}" ]] || \
        die "分类器不存在或不可读：${CLASSIFIER_QZA}"
fi

mkdir -p "${BASE_WORK}"
[[ ! -e "${WORK}" ]] || die "结果目录已存在，为防止混用旧结果而停止：${WORK}"
mkdir -p \
    "${WORK}/00_provenance" \
    "${WORK}/01_import" \
    "${WORK}/02_dada2" \
    "${WORK}/03_qc" \
    "${WORK}/04_taxonomy" \
    "${WORK}/05_phylogeny" \
    "${WORK}/06_export/asv_biom" \
    "${WORK}/06_export/rep_seqs" \
    "${WORK}/06_export/denoising_stats" \
    "${WORK}/06_export/base_transition_stats" \
    "${WORK}/06_export/asv_frequencies" \
    "${WORK}/06_export/sample_frequencies" \
    "${WORK}/logs" \
    "${WORK}/tmp"

LOG_FILE="${WORK}/logs/run.log"
exec > >(tee -a "${LOG_FILE}") 2>&1
trap 'on_error "$?" "${LINENO}" "${BASH_COMMAND}"' ERR
trap on_exit EXIT

printf '[开始] %s\n' "$(date '+%F %T')"
printf '[工作目录] %s\n' "${WORK}"

activate_qiime_environment

for executable in conda qiime biom python gzip find sort awk tail comm sha256sum xargs; do
    command -v "${executable}" >/dev/null 2>&1 || \
        die "当前环境缺少命令：${executable}。请在运行前安装，不要在正式流程中临时修改环境。"
done

QIIME_VERSION="$(qiime --version 2>&1)"
printf '%s\n' "${QIIME_VERSION}" | tee "${WORK}/00_provenance/qiime-version.txt"
if [[ -n "${EXPECTED_QIIME_VERSION}" && "${QIIME_VERSION}" != *"${EXPECTED_QIIME_VERSION}"* ]]; then
    die "当前 QIIME 2 版本不含 ${EXPECTED_QIIME_VERSION}；本脚本按 2026.4 接口编写。"
fi

biom --version > "${WORK}/00_provenance/biom-version.txt"
conda list --explicit > "${WORK}/00_provenance/conda-explicit.txt"
(qiime info > "${WORK}/00_provenance/qiime-info.txt" 2>&1) || \
    warn "qiime info 执行失败，已跳过；不影响主流程。"
cp -- "${SCRIPT_SOURCE}" "${WORK}/00_provenance/run-script.sh"

if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
    qiime tools validate "${CLASSIFIER_QZA}"
    qiime tools peek "${CLASSIFIER_QZA}" | tee "${WORK}/00_provenance/classifier-peek.txt"
    grep -q 'TaxonomicClassifier' "${WORK}/00_provenance/classifier-peek.txt" || \
        die "CLASSIFIER_QZA 不是 TaxonomicClassifier artifact。"
    sha256sum "${CLASSIFIER_QZA}" > "${WORK}/00_provenance/classifier.sha256"
    python -c 'import sklearn; print(sklearn.__version__)' \
        > "${WORK}/00_provenance/runtime-sklearn-version.txt"
fi

export TMPDIR="${WORK}/tmp"
cd "${WORK}"

R1_LIST="00_provenance/r1-files.nul"
R2_LIST="00_provenance/r2-files.nul"

find -L "${CLEAN_DIR}" -mindepth 1 -maxdepth "${MAX_INPUT_DEPTH}" \
    -type f -name '*_R1.fastq.gz' ! -name '*_raw_R1.fastq.gz' -print0 \
    | sort -z > "${R1_LIST}"
find -L "${CLEAN_DIR}" -mindepth 1 -maxdepth "${MAX_INPUT_DEPTH}" \
    -type f -name '*_R2.fastq.gz' ! -name '*_raw_R2.fastq.gz' -print0 \
    | sort -z > "${R2_LIST}"

mapfile -d '' -t R1_FILES < "${R1_LIST}"
mapfile -d '' -t R2_FILES < "${R2_LIST}"

[[ "${#R1_FILES[@]}" -gt 0 ]] || \
    die "未找到 *_R1.fastq.gz；请检查 CLEAN_DIR 和 MAX_INPUT_DEPTH。"

printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n' > manifest.tsv
declare -A SEEN_SAMPLE_IDS=()

for R1 in "${R1_FILES[@]}"; do
    R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
    [[ -f "${R2}" ]] || die "找不到与 R1 配对的 R2：${R2}"
    [[ -s "${R1}" && -r "${R1}" ]] || die "R1 为空或不可读：${R1}"
    [[ -s "${R2}" && -r "${R2}" ]] || die "R2 为空或不可读：${R2}"

    SAMPLE="$(basename "${R1%_R1.fastq.gz}")"
    [[ -n "${SAMPLE}" ]] || die "发现空 sample ID：${R1}"
    [[ "${SAMPLE}" != \#* ]] || die "sample ID 不能以 # 开头：${SAMPLE}"
    [[ "${SAMPLE}" != *$'\t'* && "${SAMPLE}" != *$'\r'* && "${SAMPLE}" != *$'\n'* ]] || \
        die "sample ID 含 tab/换行符：${SAMPLE}"
    [[ -z "${SEEN_SAMPLE_IDS[${SAMPLE}]+x}" ]] || die "发现重复 sample ID：${SAMPLE}"
    SEEN_SAMPLE_IDS["${SAMPLE}"]=1

    printf '%s\t%s\t%s\n' "${SAMPLE}" "${R1}" "${R2}" >> manifest.tsv
done

for R2 in "${R2_FILES[@]}"; do
    R1="${R2%_R2.fastq.gz}_R1.fastq.gz"
    [[ -f "${R1}" ]] || die "发现没有对应 R1 的孤立 R2：${R2}"
done

N_SAMPLES="${#R1_FILES[@]}"
[[ "${#R2_FILES[@]}" -eq "${N_SAMPLES}" ]] || \
    die "R1/R2 文件数不一致：R1=${N_SAMPLES}, R2=${#R2_FILES[@]}"

if [[ "${CHECK_GZIP}" -eq 1 ]]; then
    printf '[输入检查] 正在验证 %s 对 gzip 文件完整性...\n' "${N_SAMPLES}"
    for R1 in "${R1_FILES[@]}"; do
        R2="${R1%_R1.fastq.gz}_R2.fastq.gz"
        gzip -t -- "${R1}" "${R2}"
    done
fi

cut -f1 manifest.tsv | tail -n +2 | LC_ALL=C sort > 00_provenance/manifest-sample-ids.txt

METADATA_ARGS=()
if [[ -n "${METADATA_FILE}" ]]; then
    METADATA_ARGS=(--m-metadata-file "${METADATA_FILE}")
    awk -F '\t' 'NR > 1 && $1 !~ /^#/ && $1 != "" {sub(/\r$/, "", $1); print $1}' \
        "${METADATA_FILE}" | LC_ALL=C sort > 00_provenance/metadata-sample-ids.txt

    METADATA_DUPLICATES="$(uniq -d 00_provenance/metadata-sample-ids.txt || true)"
    [[ -z "${METADATA_DUPLICATES}" ]] || \
        die "metadata 中存在重复 sample ID：${METADATA_DUPLICATES}"

    MISSING_METADATA_IDS="$(LC_ALL=C comm -23 \
        00_provenance/manifest-sample-ids.txt \
        00_provenance/metadata-sample-ids.txt || true)"
    [[ -z "${MISSING_METADATA_IDS}" ]] || \
        die "下列 manifest 样本缺少 metadata：${MISSING_METADATA_IDS}"

    EXTRA_METADATA_IDS="$(LC_ALL=C comm -13 \
        00_provenance/manifest-sample-ids.txt \
        00_provenance/metadata-sample-ids.txt || true)"
    if [[ -n "${EXTRA_METADATA_IDS}" ]]; then
        warn "metadata 含本次 manifest 之外的样本，将被下游忽略：${EXTRA_METADATA_IDS}"
    fi

    qiime metadata tabulate \
        --m-input-file "${METADATA_FILE}" \
        --o-visualization 00_provenance/sample-metadata.qzv
fi

printf '识别到样本数：%s\n' "${N_SAMPLES}"
printf 'manifest 前 5 个样本：\n'
head -n 6 manifest.tsv

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.tsv \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path 01_import/demux.qza

qiime demux summarize \
    --i-data 01_import/demux.qza \
    --o-visualization 01_import/demux.qzv

if [[ "${TRUNC_F}" -eq 0 || "${TRUNC_R}" -eq 0 ]]; then
    warn "TRUNC_F/TRUNC_R 至少一个为 0。该设置合法，但请依据 demux.qzv 和扩增子长度检查合并保留率。"
fi

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs 01_import/demux.qza \
    --p-trim-left-f "${TRIM_LEFT_F}" \
    --p-trim-left-r "${TRIM_LEFT_R}" \
    --p-trunc-len-f "${TRUNC_F}" \
    --p-trunc-len-r "${TRUNC_R}" \
    --p-max-ee-f "${MAX_EE_F}" \
    --p-max-ee-r "${MAX_EE_R}" \
    --p-trunc-q "${TRUNC_Q}" \
    --p-min-overlap "${MIN_OVERLAP}" \
    --p-max-merge-mismatch "${MAX_MERGE_MISMATCH}" \
    --p-pooling-method "${POOLING_METHOD}" \
    --p-n-reads-learn "${N_READS_LEARN}" \
    --p-chimera-method consensus \
    --p-n-threads "${THREADS}" \
    --o-table 02_dada2/asv-table.qza \
    --o-representative-sequences 02_dada2/rep-seqs.qza \
    --o-denoising-stats 02_dada2/denoising-stats.qza \
    --o-base-transition-stats 02_dada2/base-transition-stats.qza

qiime metadata tabulate \
    --m-input-file 02_dada2/denoising-stats.qza \
    --o-visualization 03_qc/denoising-stats.qzv

qiime dada2 plot-base-transitions \
    --i-base-transition-stats 02_dada2/base-transition-stats.qza \
    --o-visualization 03_qc/base-transition-stats.qzv

qiime feature-table summarize \
    --i-table 02_dada2/asv-table.qza \
    "${METADATA_ARGS[@]}" \
    --o-feature-frequencies 03_qc/asv-frequencies.qza \
    --o-sample-frequencies 03_qc/sample-frequencies.qza \
    --o-summary 03_qc/asv-table.qzv

qiime feature-table tabulate-seqs \
    --i-data 02_dada2/rep-seqs.qza \
    --o-visualization 03_qc/rep-seqs.qzv

if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
    qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER_QZA}" \
        --i-reads 02_dada2/rep-seqs.qza \
        --p-n-jobs "${TAX_N_JOBS}" \
        --p-confidence "${TAX_CONFIDENCE}" \
        --p-read-orientation "${TAX_READ_ORIENTATION}" \
        --o-classification 04_taxonomy/taxonomy.qza

    qiime metadata tabulate \
        --m-input-file 04_taxonomy/taxonomy.qza \
        --o-visualization 04_taxonomy/taxonomy.qzv

    qiime feature-table tabulate-seqs \
        --i-data 02_dada2/rep-seqs.qza \
        --i-taxonomy 04_taxonomy/taxonomy.qza \
        --o-visualization 04_taxonomy/rep-seqs-with-taxonomy.qzv

    qiime taxa barplot \
        --i-table 02_dada2/asv-table.qza \
        --i-taxonomy 04_taxonomy/taxonomy.qza \
        "${METADATA_ARGS[@]}" \
        --o-visualization 04_taxonomy/taxa-barplot.qzv
fi

if [[ "${RUN_TREE}" -eq 1 ]]; then
    qiime alignment mafft \
        --i-sequences 02_dada2/rep-seqs.qza \
        --p-n-threads "${MAFFT_THREADS}" \
        --o-alignment 05_phylogeny/aligned-rep-seqs.qza

    qiime alignment mask \
        --i-alignment 05_phylogeny/aligned-rep-seqs.qza \
        --p-max-gap-frequency "${MASK_MAX_GAP_FREQUENCY}" \
        --p-min-conservation "${MASK_MIN_CONSERVATION}" \
        --o-masked-alignment 05_phylogeny/masked-aligned-rep-seqs.qza

    qiime phylogeny fasttree \
        --i-alignment 05_phylogeny/masked-aligned-rep-seqs.qza \
        --p-n-threads "${FASTTREE_THREADS}" \
        --o-tree 05_phylogeny/unrooted-tree.qza

    qiime phylogeny midpoint-root \
        --i-tree 05_phylogeny/unrooted-tree.qza \
        --o-rooted-tree 05_phylogeny/rooted-tree.qza
fi

qiime tools export \
    --input-path 02_dada2/asv-table.qza \
    --output-path 06_export/asv_biom

biom convert \
    -i 06_export/asv_biom/feature-table.biom \
    -o 06_export/ASV_table.biom.tsv \
    --to-tsv

awk -F '\t' 'BEGIN {OFS="\t"} /^# Constructed/ {next} $1=="#OTU ID" {$1="Feature ID"} {print}' \
    06_export/ASV_table.biom.tsv > 06_export/ASV_table.tsv

qiime tools export \
    --input-path 02_dada2/rep-seqs.qza \
    --output-path 06_export/rep_seqs
qiime tools export \
    --input-path 02_dada2/denoising-stats.qza \
    --output-path 06_export/denoising_stats
qiime tools export \
    --input-path 02_dada2/base-transition-stats.qza \
    --output-path 06_export/base_transition_stats
qiime tools export \
    --input-path 03_qc/asv-frequencies.qza \
    --output-path 06_export/asv_frequencies
qiime tools export \
    --input-path 03_qc/sample-frequencies.qza \
    --output-path 06_export/sample_frequencies

if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
    mkdir -p 06_export/taxonomy
    qiime tools export \
        --input-path 04_taxonomy/taxonomy.qza \
        --output-path 06_export/taxonomy

    {
        printf '#OTUID\ttaxonomy\tconfidence\n'
        tail -n +2 06_export/taxonomy/taxonomy.tsv
    } > 06_export/taxonomy/biom-taxonomy.tsv

    biom add-metadata \
        -i 06_export/asv_biom/feature-table.biom \
        -o 06_export/asv_biom/feature-table-with-taxonomy.biom \
        --observation-metadata-fp 06_export/taxonomy/biom-taxonomy.tsv \
        --sc-separated taxonomy

    biom convert \
        -i 06_export/asv_biom/feature-table-with-taxonomy.biom \
        -o 06_export/ASV_table_with_taxonomy.biom.tsv \
        --to-tsv \
        --header-key taxonomy

    awk -F '\t' 'BEGIN {OFS="\t"} /^# Constructed/ {next} $1=="#OTU ID" {$1="Feature ID"} {print}' \
        06_export/ASV_table_with_taxonomy.biom.tsv > 06_export/ASV_table_with_taxonomy.tsv
fi

if [[ "${RUN_TREE}" -eq 1 ]]; then
    mkdir -p \
        06_export/phylogeny/alignment \
        06_export/phylogeny/masked_alignment \
        06_export/phylogeny/unrooted_tree \
        06_export/phylogeny/rooted_tree

    qiime tools export \
        --input-path 05_phylogeny/aligned-rep-seqs.qza \
        --output-path 06_export/phylogeny/alignment
    qiime tools export \
        --input-path 05_phylogeny/masked-aligned-rep-seqs.qza \
        --output-path 06_export/phylogeny/masked_alignment
    qiime tools export \
        --input-path 05_phylogeny/unrooted-tree.qza \
        --output-path 06_export/phylogeny/unrooted_tree
    qiime tools export \
        --input-path 05_phylogeny/rooted-tree.qza \
        --output-path 06_export/phylogeny/rooted_tree
fi

cat > 00_provenance/run_parameters.txt <<PARAMS
QIIME2 environment: ${ENV_NAME}
QIIME2 version: ${QIIME_VERSION}
Input CleanData directory: ${CLEAN_DIR}
Metadata file: ${METADATA_FILE:-not provided}
Number of samples: ${N_SAMPLES}
Threads: ${THREADS}
trim-left-f: ${TRIM_LEFT_F}
trim-left-r: ${TRIM_LEFT_R}
trunc-len-f: ${TRUNC_F}
trunc-len-r: ${TRUNC_R}
max-ee-f: ${MAX_EE_F}
max-ee-r: ${MAX_EE_R}
trunc-q: ${TRUNC_Q}
min-overlap: ${MIN_OVERLAP}
max-merge-mismatch: ${MAX_MERGE_MISMATCH}
pooling-method: ${POOLING_METHOD}
n-reads-learn: ${N_READS_LEARN}
chimera-method: consensus
run-taxonomy: ${RUN_TAXONOMY}
amplicon-region: ${AMPLICON_REGION}
taxonomy-db: ${TAXONOMY_DB}
classifier: ${CLASSIFIER_QZA:-not provided}
taxonomy-n-jobs: ${TAX_N_JOBS}
taxonomy-confidence: ${TAX_CONFIDENCE}
taxonomy-read-orientation: ${TAX_READ_ORIENTATION}
run-tree: ${RUN_TREE}
mafft-threads: ${MAFFT_THREADS}
fasttree-threads: ${FASTTREE_THREADS}
mask-max-gap-frequency: ${MASK_MAX_GAP_FREQUENCY}
mask-min-conservation: ${MASK_MIN_CONSERVATION}
PARAMS

cp 00_provenance/run_parameters.txt 06_export/run_parameters.txt
find 06_export -type f -print0 | sort -z | xargs -0 sha256sum \
    > 00_provenance/exported-files.sha256
touch _SUCCESS

printf 'ASV 流程完成：%s\n' "${WORK}"
printf 'ASV 丰度表：%s/06_export/ASV_table.tsv\n' "${WORK}"
printf 'ASV 代表序列：%s/06_export/rep_seqs/dna-sequences.fasta\n' "${WORK}"
printf 'DADA2 统计：%s/06_export/denoising_stats/stats.tsv\n' "${WORK}"
if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
    printf '分类注释：%s/06_export/taxonomy/taxonomy.tsv\n' "${WORK}"
    printf '带注释丰度表：%s/06_export/ASV_table_with_taxonomy.tsv\n' "${WORK}"
fi
if [[ "${RUN_TREE}" -eq 1 ]]; then
    printf '有根树：%s/06_export/phylogeny/rooted_tree/tree.nwk\n' "${WORK}"
fi
printf '日志：%s\n' "${LOG_FILE}"
