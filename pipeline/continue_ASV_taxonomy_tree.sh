#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

ENV_NAME="${ENV_NAME:-rachis-qiime2-2026.4}"

WORK="${WORK:?Set WORK to an existing QIIME2 run directory}"

SILVA_DIR="${SILVA_DIR:-data/reference/silva}"

REF_SEQS="${REF_SEQS:-${SILVA_DIR}/silva-138.2-ssu-nr99-seqs-V3-V4-uniq.qza}"
REF_TAX="${REF_TAX:-${SILVA_DIR}/silva-138.2-ssu-nr99-tax-V3-V4-derep-uniq.qza}"

TAX_N_JOBS="${TAX_N_JOBS:-4}"
MAFFT_THREADS="${MAFFT_THREADS:-16}"
FASTTREE_THREADS="${FASTTREE_THREADS:-1}"

TAX_CONFIDENCE="${TAX_CONFIDENCE:-0.7}"
TAX_READ_ORIENTATION="${TAX_READ_ORIENTATION:-auto}"

MASK_MAX_GAP_FREQUENCY="${MASK_MAX_GAP_FREQUENCY:-1.0}"
MASK_MIN_CONSERVATION="${MASK_MIN_CONSERVATION:-0.4}"

START_DIR="${PWD}"
for path_name in WORK REF_SEQS REF_TAX; do
    path_value="${!path_name}"
    if [[ "${path_value}" != /* ]]; then
        printf -v "${path_name}" '%s/%s' "${START_DIR}" "${path_value}"
    fi
done
RESULT_DIR="${RESULT_DIR:-${WORK}/taxonomy_tree}"
if [[ "${RESULT_DIR}" != /* ]]; then
    RESULT_DIR="${START_DIR}/${RESULT_DIR}"
fi
if [[ -e "${RESULT_DIR}" ]]; then
    printf 'Output path already exists: %s\nChoose a new RESULT_DIR.\n' "${RESULT_DIR}" >&2
    exit 1
fi
CLASSIFIER_DIR="${CLASSIFIER_DIR:-${RESULT_DIR}/reference_classifier}"
if [[ "${CLASSIFIER_DIR}" != /* ]]; then
    CLASSIFIER_DIR="${START_DIR}/${CLASSIFIER_DIR}"
fi

source "${CONDA_SH:-$(conda info --base)/etc/profile.d/conda.sh}"
conda activate "${ENV_NAME}"

cd "${WORK}" || {
    echo "[错误] 找不到结果目录：${WORK}"
    exit 1
}

echo "当前工作目录：${WORK}"
echo "QIIME2 版本："
qiime --version

SKLEARN_VERSION="$(python -c 'import sklearn; print(sklearn.__version__)')"
echo "当前 scikit-learn 版本：${SKLEARN_VERSION}"

ASV_TABLE="${WORK}/02_dada2/asv-table.qza"
REP_SEQS="${WORK}/02_dada2/rep-seqs.qza"
DENOISING_STATS="${WORK}/02_dada2/denoising-stats.qza"
BASE_TRANSITION_STATS="${WORK}/02_dada2/base-transition-stats.qza"

for file in \
    "${ASV_TABLE}" \
    "${REP_SEQS}" \
    "${DENOISING_STATS}" \
    "${BASE_TRANSITION_STATS}" \
    "${WORK}/03_qc/asv-frequencies.qza" \
    "${WORK}/03_qc/sample-frequencies.qza"; do

    if [[ ! -s "${file}" ]]; then
        echo "[错误] 文件不存在或为空：${file}"
        exit 1
    fi
done

echo "检查已有 QIIME2 ASV 文件："

qiime tools validate "${ASV_TABLE}"
qiime tools validate "${REP_SEQS}"
qiime tools validate "${DENOISING_STATS}"
qiime tools validate "${BASE_TRANSITION_STATS}"

for file in "${REF_SEQS}" "${REF_TAX}"; do
    if [[ ! -s "${file}" ]]; then
        echo "[错误] SILVA 文件不存在或为空：${file}"
        exit 1
    fi
done

echo "检查 SILVA V3-V4 参考序列："
qiime tools validate "${REF_SEQS}"
qiime tools peek "${REF_SEQS}"

echo "检查 SILVA V3-V4 分类表："
qiime tools validate "${REF_TAX}"
qiime tools peek "${REF_TAX}"

mkdir -p "${CLASSIFIER_DIR}"
mkdir -p "${RESULT_DIR}/00_provenance"

NEW_CLASSIFIER="${CLASSIFIER_DIR}/silva-138.2-ssu-nr99-V3-V4-sklearn-${SKLEARN_VERSION}.qza"

if [[ -f "${NEW_CLASSIFIER}" ]]; then
    echo "发现已有分类器，正在验证：${NEW_CLASSIFIER}"

    if qiime tools validate "${NEW_CLASSIFIER}"; then
        echo "已有分类器验证通过，直接使用。"
    else
        printf 'Existing classifier failed validation: %s\nChoose a new CLASSIFIER_DIR.\n' "${NEW_CLASSIFIER}" >&2
        exit 1
    fi
fi

if [[ ! -f "${NEW_CLASSIFIER}" ]]; then
    echo "开始训练 SILVA 138.2 V3-V4 分类器"
    echo "输出：${NEW_CLASSIFIER}"

    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads "${REF_SEQS}" \
        --i-reference-taxonomy "${REF_TAX}" \
        --o-classifier "${NEW_CLASSIFIER}"
fi

qiime tools validate "${NEW_CLASSIFIER}"
qiime tools peek "${NEW_CLASSIFIER}" \
    | tee "${RESULT_DIR}/00_provenance/new-classifier-peek.txt"

sha256sum "${NEW_CLASSIFIER}" \
    > "${RESULT_DIR}/00_provenance/new-classifier.sha256"

printf '%s\n' "${SKLEARN_VERSION}" \
    > "${RESULT_DIR}/00_provenance/classifier-sklearn-version.txt"

printf '%s\n' "${NEW_CLASSIFIER}" \
    > "${RESULT_DIR}/00_provenance/classifier-path.txt"

mkdir -p "${RESULT_DIR}/04_taxonomy"

echo "开始 ASV 分类注释"

qiime feature-classifier classify-sklearn \
    --i-classifier "${NEW_CLASSIFIER}" \
    --i-reads "${REP_SEQS}" \
    --p-n-jobs "${TAX_N_JOBS}" \
    --p-confidence "${TAX_CONFIDENCE}" \
    --p-read-orientation "${TAX_READ_ORIENTATION}" \
    --o-classification "${RESULT_DIR}/04_taxonomy/taxonomy.qza"

qiime metadata tabulate \
    --m-input-file "${RESULT_DIR}/04_taxonomy/taxonomy.qza" \
    --o-visualization "${RESULT_DIR}/04_taxonomy/taxonomy.qzv"

qiime feature-table tabulate-seqs \
    --i-data "${REP_SEQS}" \
    --i-taxonomy "${RESULT_DIR}/04_taxonomy/taxonomy.qza" \
    --o-visualization "${RESULT_DIR}/04_taxonomy/rep-seqs-with-taxonomy.qzv"

qiime taxa barplot \
    --i-table "${ASV_TABLE}" \
    --i-taxonomy "${RESULT_DIR}/04_taxonomy/taxonomy.qza" \
    --o-visualization "${RESULT_DIR}/04_taxonomy/taxa-barplot.qzv"

mkdir -p "${RESULT_DIR}/05_phylogeny"

echo "开始构建系统发育树"

qiime alignment mafft \
    --i-sequences "${REP_SEQS}" \
    --p-n-threads "${MAFFT_THREADS}" \
    --o-alignment "${RESULT_DIR}/05_phylogeny/aligned-rep-seqs.qza"

qiime alignment mask \
    --i-alignment "${RESULT_DIR}/05_phylogeny/aligned-rep-seqs.qza" \
    --p-max-gap-frequency "${MASK_MAX_GAP_FREQUENCY}" \
    --p-min-conservation "${MASK_MIN_CONSERVATION}" \
    --o-masked-alignment "${RESULT_DIR}/05_phylogeny/masked-aligned-rep-seqs.qza"

qiime phylogeny fasttree \
    --i-alignment "${RESULT_DIR}/05_phylogeny/masked-aligned-rep-seqs.qza" \
    --p-n-threads "${FASTTREE_THREADS}" \
    --o-tree "${RESULT_DIR}/05_phylogeny/unrooted-tree.qza"

qiime phylogeny midpoint-root \
    --i-tree "${RESULT_DIR}/05_phylogeny/unrooted-tree.qza" \
    --o-rooted-tree "${RESULT_DIR}/05_phylogeny/rooted-tree.qza"

echo "开始导出最终结果"

mkdir -p "${RESULT_DIR}/06_export"

mkdir -p \
    "${RESULT_DIR}/06_export/asv_biom" \
    "${RESULT_DIR}/06_export/rep_seqs" \
    "${RESULT_DIR}/06_export/denoising_stats" \
    "${RESULT_DIR}/06_export/base_transition_stats" \
    "${RESULT_DIR}/06_export/asv_frequencies" \
    "${RESULT_DIR}/06_export/sample_frequencies" \
    "${RESULT_DIR}/06_export/taxonomy" \
    "${RESULT_DIR}/06_export/phylogeny/alignment" \
    "${RESULT_DIR}/06_export/phylogeny/masked_alignment" \
    "${RESULT_DIR}/06_export/phylogeny/unrooted_tree" \
    "${RESULT_DIR}/06_export/phylogeny/rooted_tree"

qiime tools export \
    --input-path "${ASV_TABLE}" \
    --output-path "${RESULT_DIR}/06_export/asv_biom"

biom convert \
    -i "${RESULT_DIR}/06_export/asv_biom/feature-table.biom" \
    -o "${RESULT_DIR}/06_export/ASV_table.biom.tsv" \
    --to-tsv

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
/^# Constructed/ {
    next
}
$1 == "#OTU ID" {
    $1 = "Feature ID"
}
{
    print
}
' "${RESULT_DIR}/06_export/ASV_table.biom.tsv" \
  > "${RESULT_DIR}/06_export/ASV_table.tsv"

qiime tools export \
    --input-path "${REP_SEQS}" \
    --output-path "${RESULT_DIR}/06_export/rep_seqs"

qiime tools export \
    --input-path "${DENOISING_STATS}" \
    --output-path "${RESULT_DIR}/06_export/denoising_stats"

qiime tools export \
    --input-path "${BASE_TRANSITION_STATS}" \
    --output-path "${RESULT_DIR}/06_export/base_transition_stats"

qiime tools export \
    --input-path "${WORK}/03_qc/asv-frequencies.qza" \
    --output-path "${RESULT_DIR}/06_export/asv_frequencies"

qiime tools export \
    --input-path "${WORK}/03_qc/sample-frequencies.qza" \
    --output-path "${RESULT_DIR}/06_export/sample_frequencies"

qiime tools export \
    --input-path "${RESULT_DIR}/04_taxonomy/taxonomy.qza" \
    --output-path "${RESULT_DIR}/06_export/taxonomy"

{
    printf '#OTUID\ttaxonomy\tconfidence\n'
    tail -n +2 "${RESULT_DIR}/06_export/taxonomy/taxonomy.tsv"
} > "${RESULT_DIR}/06_export/taxonomy/biom-taxonomy.tsv"

biom add-metadata \
    -i "${RESULT_DIR}/06_export/asv_biom/feature-table.biom" \
    -o "${RESULT_DIR}/06_export/asv_biom/feature-table-with-taxonomy.biom" \
    --observation-metadata-fp "${RESULT_DIR}/06_export/taxonomy/biom-taxonomy.tsv" \
    --sc-separated taxonomy

biom convert \
    -i "${RESULT_DIR}/06_export/asv_biom/feature-table-with-taxonomy.biom" \
    -o "${RESULT_DIR}/06_export/ASV_table_with_taxonomy.biom.tsv" \
    --to-tsv \
    --header-key taxonomy

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
/^# Constructed/ {
    next
}
$1 == "#OTU ID" {
    $1 = "Feature ID"
}
{
    print
}
' "${RESULT_DIR}/06_export/ASV_table_with_taxonomy.biom.tsv" \
  > "${RESULT_DIR}/06_export/ASV_table_with_taxonomy.tsv"

qiime tools export \
    --input-path "${RESULT_DIR}/05_phylogeny/aligned-rep-seqs.qza" \
    --output-path "${RESULT_DIR}/06_export/phylogeny/alignment"

qiime tools export \
    --input-path "${RESULT_DIR}/05_phylogeny/masked-aligned-rep-seqs.qza" \
    --output-path "${RESULT_DIR}/06_export/phylogeny/masked_alignment"

qiime tools export \
    --input-path "${RESULT_DIR}/05_phylogeny/unrooted-tree.qza" \
    --output-path "${RESULT_DIR}/06_export/phylogeny/unrooted_tree"

qiime tools export \
    --input-path "${RESULT_DIR}/05_phylogeny/rooted-tree.qza" \
    --output-path "${RESULT_DIR}/06_export/phylogeny/rooted_tree"

cat > "${RESULT_DIR}/06_export/taxonomy_run_parameters.txt" <<PARAMS
QIIME2 environment: ${ENV_NAME}
Current scikit-learn version: ${SKLEARN_VERSION}
SILVA database: SILVA 138.2 SSU NR99
Amplicon region: V3-V4
Reference sequences: ${REF_SEQS}
Reference taxonomy: ${REF_TAX}
Classifier: ${NEW_CLASSIFIER}
Taxonomy confidence: ${TAX_CONFIDENCE}
Taxonomy read orientation: ${TAX_READ_ORIENTATION}
Taxonomy jobs: ${TAX_N_JOBS}
MAFFT threads: ${MAFFT_THREADS}
FastTree threads: ${FASTTREE_THREADS}
Mask max gap frequency: ${MASK_MAX_GAP_FREQUENCY}
Mask min conservation: ${MASK_MIN_CONSERVATION}
PARAMS

find "${RESULT_DIR}/06_export" -type f -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${RESULT_DIR}/00_provenance/final-exported-files.sha256"

touch "${RESULT_DIR}/_SUCCESS"

echo
echo "全部处理完成"
echo "结果目录：${RESULT_DIR}"
echo
echo "ASV 丰度表："
echo "${RESULT_DIR}/06_export/ASV_table.tsv"
echo
echo "带分类注释的 ASV 丰度表："
echo "${RESULT_DIR}/06_export/ASV_table_with_taxonomy.tsv"
echo
echo "ASV 代表序列："
echo "${RESULT_DIR}/06_export/rep_seqs/dna-sequences.fasta"
echo
echo "DADA2 统计表："
echo "${RESULT_DIR}/06_export/denoising_stats/stats.tsv"
echo
echo "分类注释表："
echo "${RESULT_DIR}/06_export/taxonomy/taxonomy.tsv"
echo
echo "物种分类可视化："
echo "${RESULT_DIR}/04_taxonomy/taxonomy.qzv"
echo
echo "物种组成柱状图："
echo "${RESULT_DIR}/04_taxonomy/taxa-barplot.qzv"
echo
echo "有根系统发育树："
echo "${RESULT_DIR}/06_export/phylogeny/rooted_tree/tree.nwk"
echo
echo "新训练分类器："
echo "${NEW_CLASSIFIER}"
