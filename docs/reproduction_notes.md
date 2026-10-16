# Reproduction notes

These notes record unresolved provenance and input requirements. They should be resolved against the original run records before a public, DOI-tagged release.

## Single rarefaction and candidate subset

`01_分析项目/00_数据预处理与输入审计/分析代码/运行代码/抽平并分发ASV表.R`

- Expects the historical runtime layout (tpm.csv and 致病菌/tpm.csv); distribution overwrites runtime copies. Do not point it at source data folders.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Taxonomy conversion

`01_分析项目/01_分类注释_SILVA_Greengenes/分析代码/运行代码/整理双注释表.R`

- NC Methods uses SILVA 138.2. Greengenes is historical comparison, not the selected primary annotation.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Alpha diversity calculation

`01_分析项目/03_Alpha多样性/分析代码/运行代码/整体_Alpha计算.R`

- Calls rrarefy at the minimum input depth. For the manuscript the input must already be the locked 21000-read complete-community table; equal depths make this step a no-op.
- Shannon is calculated with base=2; retain this statistical definition.

## Fig 2C richness GAMM calculation

`01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/整体_Alpha年际趋势.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 6A candidate time trends

`01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/致病菌_Richness与相对丰度总占比_水位子图.R`

- Use community_rarefied.csv and candidate row-subset pathogen_rarefied.csv; never use the independently rarefied candidate alpha1.csv.

## Supplementary Fig 6 site PCoA

`01_分析项目/05_PCoA分组/分析代码/运行代码/整体_PCoA分组.R`

- Run using whole-community rarefied input and explicit group tables.

## Supplementary Fig 13 candidate PCoA

`01_分析项目/05_PCoA分组/分析代码/运行代码/致病菌_PCoA分组.R`

- Candidate input must be recorded explicitly; old runtime files may differ from the once-rarefied complete-community subset.

## PCoA1 time course and RF response preparation

`01_分析项目/06_PCoA1年际变化/分析代码/运行代码/整体_PCoA1年际变化.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## PCoA GAMM helper bridge

`01_分析项目/06_PCoA1年际变化/分析代码/运行代码/GAMM曲线工具.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 2A complete-community Venn

`01_分析项目/07_Venn图/分析代码/运行代码/整体_Venn.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 12 candidate Venn

`01_分析项目/07_Venn图/分析代码/运行代码/致病菌_Venn.R`

- Relative tpm.csv must be the original 2755-ASV count table, not the rarefied candidate subset.

## Fig 2A radial taxonomy network

`01_分析项目/09_物种组成/分析代码/运行代码/分类网络.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 11 candidate relative abundance

`01_分析项目/09_物种组成/分析代码/运行代码/候选致病菌相对丰度_大小粒级_单图.R`

- Verify denominators remain the complete-community rarefied totals; this is not the within-candidate genus normalization.

## MBPD grouped composition functions

`01_分析项目/09_物种组成/分析代码/运行代码/致病菌_MBPD七层级分组组成.R`

- This entry was changed on Sep-18 to use candidate subset after rarefaction. Its default outputs do not reproduce the legacy raw-count genus figure retained in the NC manuscript. Use the explicit raw-input export wrapper for Data table 23.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 2B and Supplementary Figs 3-4 monthly composition

`01_分析项目/10_物种组成年际变化/分析代码/运行代码/所有物种组成_抽平逐月堆叠柱图.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 4F-G static networks and node removal

`01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_出现率20pct_增强版.R`

- Static robustness/retention uses unweighted natural connectivity; do not replace it with the sliding-window weighted metric.

## Fig 4G node-removal figure

`01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_ConF_USF_稳定性组合图.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Static network topology figure

`01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_ConF_USF_拓扑柱状图.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Static-network legend

`01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络图例.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 4C-E sliding-window networks

`01_分析项目/12_网络年际变化/分析代码/运行代码/整体_滑窗网络_出现率20pct_CNS.R`

- Windows advance over sampled months, not a complete calendar; the Nov-centred window includes Oct-Nov-Jan.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Network GAMM helper bridge

`01_分析项目/12_网络年际变化/分析代码/运行代码/GAMM曲线工具.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 4C-E final large-axis figures

`replot_three_network_trends_large_axes.R`

- Depends on replot_inputs copies; preserve their provenance to the 20%-prevalence sliding-network outputs.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 2 LEfSe

`01_分析项目/13_LEfSe/分析代码/运行代码/整体_LEfSe.R`

- The historical output field P.adj was not actually multiple-test adjusted in this run; retain its interpretation as documented in supplement.

## Fig 3D-F RDA calculation

`work/rda_cca_noQ/run_rda_cca_noQ.R`

- Requires work/rda_cca_noQ/input/{whole,candidate,env,groupsize}.csv. The script also calculates candidate and CCA outputs not displayed in main Figure 3.
- RDA retains taxa present in >=10% of each subset before Hellinger transform; record this in run documentation.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 3D-F final full-arrow RDA plotting

`work/rda_all_arrows_label_repel_20260902/plot_rda_all_arrows_label_repel.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 3G Mantel and environment correlations

`01_分析项目/15_Mantel检验/分析代码/运行代码/整体_Mantel.R`

- Current hash differs from Aug-28 run manifest. Source diff shows only plotting changes relative to saved work copy; exact Aug-28 statistical reproduction still requires the archived input and package environment.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 3A-C environmental/spatial/temporal decay

`01_分析项目/16_距离衰减/分析代码/运行代码/整体_距离衰减.R`

- Use the Aug-28 19-variable env.csv, retaining spatial coordinates; Q excluded by code.
- USF environment Mantel P=.0013 and time P=.0011; NC prose all P<=.001 is slightly inaccurate, do not alter computed values to fit prose.

## Supplementary Fig 10 / Data table 15 RF

`outputs/20261002_manuscript_revision/RF/rerun_RF_without_Q.R`

- Requires predecessor 17 script to exist for provenance hash even though not sourced; preserve a provenance copy or remove only that assertion in the publication adaptation.
- Uses 99_运行环境/串行运行工作区/SILVA/{env.csv,PCOA1年际变化/1.csv,随机森林/groupsize.csv}; these input data must be relocated, not the R library.
- Do not interpret negative pooled All OOB variance as predictive ability.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 10 RF plotting

`outputs/20261002_manuscript_revision/RF/plot_RF_without_Q.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## RF predecessor provenance

`01_分析项目/17_随机森林/分析代码/运行代码/整体_随机森林.R`

- Historical code can include Q as a twentieth numeric predictor; do not use it as the main current RF entry.

## Supplementary Figs 7 and 16 sample-level water regression

`01_分析项目/18_Richness水位回归/分析代码/运行代码/样本水平_水位回归_All_ConF_USF.R`

- The candidate-richness branch reads old candidate alpha1.csv and is superseded by S15_candidate_richness_single_rarefaction.R. Exclude that branch from final S15 production.

## Supplementary Fig 15 / Data table 26 candidate richness regression

`work_artifact/revision_20261002/S15_candidate_richness_single_rarefaction.R`

- Accepts overall_input_dir, pathogen_input_dir and output_dir CLI arguments; documented original relative layout locates complete tpm_rarefied.csv two parent levels above candidate input.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 9 / Data table 14 VPA

`01_分析项目/33_VPA时间环境空间/分析代码/运行代码/整体_时间环境空间_VPA.R`

- Uses 19-variable candidate environment table with VIF screening; select final dbRDA-Bray results for manuscript, not its parallel Hellinger RDA results.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Fig 9 corrected VPA labels

`work_artifact/revision_20261002/S9_VPA_ConF_USF_labels.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 4A-B iCAMP large-font exports

`work_artifact/replot_icamp_large_font.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 4A final frame and legend

`work_artifact/replot_fig4A_closed_frame_wrapped_legend.R`

- Local absolute paths and local R library references must be made portable in the publication copy.

## Supplementary Data full composition and alpha exports

`work/supplement_revision/export_complete_results.R`

- Dependencies are evaluated using parse/eval rather than source; preserve function and top-level assignment names.
- This explicit wrapper restores raw candidate abundance for genus composition, matching NC exception; direct latest grouped-composition defaults do not.
- Local absolute paths and local R library references must be made portable in the publication copy.

## Fig 6A cached-layout source

`work/pathogen_dual_axis_trend/plot_pathogen_dual_axis_trend.R`

- CLI requires input_dir and output_dir. This version is a separate hash from canonical 04/致病菌_Richness与相对丰度总占比_水位子图.R; preserve its use for downstream drawing fragment.
- Water-level inset reads precomputed water GAMM tables generated by 04/水位年际趋势.R.

## Fig 6A wider plotting helper

`work_artifact/pathogen_wider.R`

- Uses text extraction followed by eval/parse; reformatting dependency files may break grep-based boundary detection.
- Dynamic readLines/eval(parse()) source-block dependency. Placeholder p_smooth=0.001 is used for drawing but final aligned layout hardcodes <0.001; report actual cached statistics.

## Fig 6A-B final aligned panels

`work_artifact/align_pathogen_pair.R`

- Contains fixed n_analyzed=476 and n_robust=181 for this dataset; preserve as figure metadata, not general analysis defaults.
- Uses text extraction followed by eval/parse; reformatting dependency files may break grep-based boundary detection.
- Uses nested readLines/eval(parse()) extraction with regex anchors in three upstream scripts; generic format conversion can break fragile anchors.

## Figure 1A: Baihe discharge source preparation

`outputs/baihe_monthly_indices_20260921/prepare_and_verify.py`

- Recent daily inputs live outside DJK16s; include distributable input copies or document source access.
- CLI has prepare/verify modes; verification pertains to preliminary historical-calibration outputs.

## Figure 1A: complete monthly grid preparation

`outputs/baihe_monthly_indices_20260921/calculate_indices.R`

- This intermediate script also computes historical-only SRI and endpoint-reversal events that are NOT the published final indices. Final figures must use Gamma_MLE from the matrix script.

## Figure 1A: calendar-month Gamma MLE SRI-1 and SDFAI

`outputs/baihe_sri_matrix_20260922/calculate.R`

- Also computes ten alternative distribution/index methods; only Gamma_MLE and SDFAI are the reported principal analysis.
- Fixed output folder is assumed to exist; remove local library path and parameterize root in publication copy.

## Figure 1A: final SDFAI export

`outputs/baihe_gamma_sdfai_20260922/export_plot.R`

- Uses --boxed for final frame. The red September-2025 arrow in the manuscript is a downstream manual layout annotation, not generated in this script.

## Figure 1C-F: monthly HRT, TP, TN and chlorophyll a

`01_分析项目/环境变化曲线/plot_monthly_environmental_curves_CNS4panel.R`

- Use direct script invocation or explicit command-line input/output arguments. The drawing spline bridges missing months while the data retain NA; describe smoothing as display-only.

## Figure 1B: monthly water level and net inflow

`work/water_level_interannual_curve/plot_water_level_interannual_curve.R`

- Hydrology CSV is generated by 37/01_构建更新环境表.R; root-discovery assumes original two-level directory placement.

## Figure 1F vector export

`work_artifact/chla_unit_fix_20260930/export_chla_svg.R`

- source() enters a script that detects --file from the wrapper; direct use needs explicit input/output arguments. Prefer main environmental plot entry with arguments.

## Figure 5B,C,E: water level, Chla and cyanobacteria CCF

`01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/运行_CCF_两两关系.R`

- The manuscript uses analysis_type=level_sensitivity; do not substitute differenced outputs or 26_CCF三个月滞后.
- Deletes two obsolete plot files within its output directory.

## Figure 5D,F,G: TN:TP and temperature CCF

`01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/补充_TOC_WT_原始水平CCF.R`

- Run separately with CCF_SUPPLEMENT_MODE=TOC_WT, AIR_TEMP and N_P_RATIO. All three modes are needed; default alone misses AIRT and TNTP.
- Current figure reports unadjusted maximum-lag permutation P as explicitly captioned; preserve method rather than silently changing to BH.

## Figure 5 original-level CCF panel plotting

`01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/绘图_CCF_六关系_统一优化.R`

- Reads input/output environment variables CCF_LAG_FILE, CCF_PEAK_FILE etc.; project_dir derived from script placement.

## Figure 5A: monthly cyanobacterial abundance

`01_分析项目/34_蓝藻门结构方程模型/分析代码/03_蓝藻门月均相对丰度纵向曲线.R`

- Located inside SEM directory but this plotting script is part of final paper; SEM model-fitting scripts are not.

## Figure 5A local launcher

`run_cyanobacteria_vertical_curve.R`

- Local R library path only; redundant once a portable entrypoint is supplied.

## Figure 5 final vector-panel layout

`work_artifact/layout_figure5_aligned.py`

- Requires publication panel PDFs under their copied final names; no data refitting.

## Figure 6B,C: ASV-level TITAN2

`01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_致病菌候选ASV_抽平数据.R`

- Use WL formal; smoke mode is not final. Local absolute project root and R-library path must be parameterized.

## Figure 6C: top 20 unique-species representative ASVs

`01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_WL_top20_不同种代表稳健ASV山脊图_最新注释.R`

- Final publication label abbreviations require the 2026-09-28 helper with --unambiguous.

## Figure 6C: final unambiguous genus labels

`work_artifact/relabel_titan_genus_20260928.R`

- Required argument --unambiguous. Uses readLines/eval(parse()) to extract drawing block from titan_ridge via anchors; cleanup must preserve anchors or replace dynamic loading safely.

## Local TITAN ridge launcher

`replot_fig6b_titan_labels.R`

- Superseded for final labels by 2026-09-28 unambiguous helper; wrapper itself is optional.

## Supplementary Fig. 11: full richness and candidate fraction

`01_分析项目/27_结果汇总/分析代码/运行代码/整体_Richness与候选致病菌相对丰度组合图.R`

- Requires five arguments: complete-community alpha CSV, group CSV, fraction-abundance CSV, figure directory, table directory. Alpha CSV should be current 21000-read table from analysis 03.

## Supplementary Fig. 1: four-station study map

`figures/study_area/plot_accurate_study_area.R`

- All shapefile sidecars must be retained. OSM-derived water geometry needs source/attribution and redistribution terms documented.
- No separate source scripts exist for _polished/_CNS filenames; do not infer unique versions from those names.
