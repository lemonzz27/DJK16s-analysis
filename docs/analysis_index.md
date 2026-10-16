# Analysis index

Paths retain the source project layout so that local dependencies remain resolvable. Use the launcher and working-directory catalogue in the run guide for R scripts; Python preparation/layout commands start from the repository root. `confirmed` means a match to manuscript methods, values or figure provenance; it does not mean the complete analysis was rerun during packaging.

| Analysis | Script | Role |
| --- | --- | --- |
| Single rarefaction and candidate subset | [抽平并分发ASV表.R](../01_分析项目/00_数据预处理与输入审计/分析代码/运行代码/抽平并分发ASV表.R) | confirmed |
| Taxonomy conversion | [整理双注释表.R](../01_分析项目/01_分类注释_SILVA_Greengenes/分析代码/运行代码/整理双注释表.R) | helper |
| Alpha diversity calculation | [整体_Alpha计算.R](../01_分析项目/03_Alpha多样性/分析代码/运行代码/整体_Alpha计算.R) | confirmed |
| Supplementary Fig 5 alpha comparison | [整体_Alpha分组箱线图.R](../01_分析项目/03_Alpha多样性/分析代码/运行代码/整体_Alpha分组箱线图.R) | confirmed |
| Shared GAMM functions | [GAMM曲线工具.R](../01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/GAMM曲线工具.R) | helper |
| Fig 2C richness GAMM calculation | [整体_Alpha年际趋势.R](../01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/整体_Alpha年际趋势.R) | confirmed |
| Water-level GAMM inset | [水位年际趋势.R](../01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/水位年际趋势.R) | confirmed |
| Fig 2C richness presentation | [整体_Richness与水位年际趋势组合.R](../01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/整体_Richness与水位年际趋势组合.R) | confirmed |
| Fig 6A candidate time trends | [致病菌_Richness与相对丰度总占比_水位子图.R](../01_分析项目/04_Alpha多样性年际变化/分析代码/运行代码/致病菌_Richness与相对丰度总占比_水位子图.R) | confirmed |
| Fig 2D month and size PCoA | [PCoA_月份大小合并.R](../01_分析项目/05_PCoA分组/分析代码/运行代码/PCoA_月份大小合并.R) | confirmed |
| Supplementary Fig 6 site PCoA | [整体_PCoA分组.R](../01_分析项目/05_PCoA分组/分析代码/运行代码/整体_PCoA分组.R) | confirmed |
| Supplementary Fig 13 candidate PCoA | [致病菌_PCoA分组.R](../01_分析项目/05_PCoA分组/分析代码/运行代码/致病菌_PCoA分组.R) | confirmed |
| PCoA1 time course and RF response preparation | [整体_PCoA1年际变化.R](../01_分析项目/06_PCoA1年际变化/分析代码/运行代码/整体_PCoA1年际变化.R) | confirmed |
| PCoA GAMM helper bridge | [GAMM曲线工具.R](../01_分析项目/06_PCoA1年际变化/分析代码/运行代码/GAMM曲线工具.R) | helper |
| Fig 2A complete-community Venn | [整体_Venn.R](../01_分析项目/07_Venn图/分析代码/运行代码/整体_Venn.R) | confirmed |
| Supplementary Fig 12 candidate Venn | [致病菌_Venn.R](../01_分析项目/07_Venn图/分析代码/运行代码/致病菌_Venn.R) | confirmed |
| Fig 2A radial taxonomy network | [分类网络.R](../01_分析项目/09_物种组成/分析代码/运行代码/分类网络.R) | confirmed |
| Supplementary Fig 11 candidate relative abundance | [候选致病菌相对丰度_大小粒级_单图.R](../01_分析项目/09_物种组成/分析代码/运行代码/候选致病菌相对丰度_大小粒级_单图.R) | confirmed |
| MBPD grouped composition functions | [致病菌_MBPD七层级分组组成.R](../01_分析项目/09_物种组成/分析代码/运行代码/致病菌_MBPD七层级分组组成.R) | helper |
| Fig 2B and Supplementary Figs 3-4 monthly composition | [所有物种组成_抽平逐月堆叠柱图.R](../01_分析项目/10_物种组成年际变化/分析代码/运行代码/所有物种组成_抽平逐月堆叠柱图.R) | confirmed |
| Fig 4F-G static networks and node removal | [静态网络_出现率20pct_增强版.R](../01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_出现率20pct_增强版.R) | confirmed |
| Fig 4G node-removal figure | [静态网络_ConF_USF_稳定性组合图.R](../01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_ConF_USF_稳定性组合图.R) | confirmed |
| Static network topology figure | [静态网络_ConF_USF_拓扑柱状图.R](../01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络_ConF_USF_拓扑柱状图.R) | helper |
| Static-network legend | [静态网络图例.R](../01_分析项目/11_静态共现网络/分析代码/运行代码/静态网络图例.R) | helper |
| Fig 4C-E sliding-window networks | [整体_滑窗网络_出现率20pct_CNS.R](../01_分析项目/12_网络年际变化/分析代码/运行代码/整体_滑窗网络_出现率20pct_CNS.R) | confirmed |
| Network GAMM helper bridge | [GAMM曲线工具.R](../01_分析项目/12_网络年际变化/分析代码/运行代码/GAMM曲线工具.R) | helper |
| Fig 4C-E final large-axis figures | [replot_three_network_trends_large_axes.R](../replot_three_network_trends_large_axes.R) | confirmed |
| Supplementary Fig 2 LEfSe | [整体_LEfSe.R](../01_分析项目/13_LEfSe/分析代码/运行代码/整体_LEfSe.R) | confirmed |
| Fig 3D-F RDA calculation | [run_rda_cca_noQ.R](../work/rda_cca_noQ/run_rda_cca_noQ.R) | confirmed |
| Fig 3D-F final full-arrow RDA plotting | [plot_rda_all_arrows_label_repel.R](../work/rda_all_arrows_label_repel_20260902/plot_rda_all_arrows_label_repel.R) | confirmed |
| Fig 3G Mantel and environment correlations | [整体_Mantel.R](../01_分析项目/15_Mantel检验/分析代码/运行代码/整体_Mantel.R) | confirmed |
| Fig 3A-C environmental/spatial/temporal decay | [整体_距离衰减.R](../01_分析项目/16_距离衰减/分析代码/运行代码/整体_距离衰减.R) | confirmed |
| Fig 3A-C combined plot | [整体_距离衰减_All_Large_Small组合.R](../01_分析项目/16_距离衰减/分析代码/运行代码/整体_距离衰减_All_Large_Small组合.R) | helper |
| Supplementary Fig 8 water-level decay | [water_level_distance_decay.R](../01_分析项目/16_距离衰减/分析代码/运行代码/water_level_distance_decay.R) | confirmed |
| Supplementary Fig 10 / Data table 15 RF | [rerun_RF_without_Q.R](../outputs/20261002_manuscript_revision/RF/rerun_RF_without_Q.R) | confirmed |
| Supplementary Fig 10 RF plotting | [plot_RF_without_Q.R](../outputs/20261002_manuscript_revision/RF/plot_RF_without_Q.R) | helper |
| RF predecessor provenance | [整体_随机森林.R](../01_分析项目/17_随机森林/分析代码/运行代码/整体_随机森林.R) | helper |
| Supplementary Figs 7 and 16 sample-level water regression | [样本水平_水位回归_All_ConF_USF.R](../01_分析项目/18_Richness水位回归/分析代码/运行代码/样本水平_水位回归_All_ConF_USF.R) | confirmed |
| Supplementary Fig 15 / Data table 26 candidate richness regression | [S15_candidate_richness_single_rarefaction.R](../work_artifact/revision_20261002/S15_candidate_richness_single_rarefaction.R) | confirmed |
| Supplementary Fig 9 / Data table 14 VPA | [整体_时间环境空间_VPA.R](../01_分析项目/33_VPA时间环境空间/分析代码/运行代码/整体_时间环境空间_VPA.R) | confirmed |
| Supplementary Fig 9 corrected VPA labels | [S9_VPA_ConF_USF_labels.R](../work_artifact/revision_20261002/S9_VPA_ConF_USF_labels.R) | helper |
| Fig 4A iCAMP aggregation | [plot_icamp_assembly.R](../01_分析项目/icamp_assembly_plots/plot_icamp_assembly.R) | confirmed |
| Fig 4B stochastic/deterministic GAMM | [plot_stochastic_deterministic_ratio_trend.R](../01_分析项目/icamp_assembly_plots/plot_stochastic_deterministic_ratio_trend.R) | confirmed |
| Fig 4A-B iCAMP large-font exports | [replot_icamp_large_font.R](../work_artifact/replot_icamp_large_font.R) | helper |
| Fig 4A final frame and legend | [replot_fig4A_closed_frame_wrapped_legend.R](../work_artifact/replot_fig4A_closed_frame_wrapped_legend.R) | helper |
| Supplementary Data full composition and alpha exports | [export_complete_results.R](../work/supplement_revision/export_complete_results.R) | confirmed |
| Original radial-network calculation helper | [plot_taxonomic_network.R](../work/taxonomic_network_palette/plot_taxonomic_network.R) | helper |
| Fig 6A cached-layout source | [plot_pathogen_dual_axis_trend.R](../work/pathogen_dual_axis_trend/plot_pathogen_dual_axis_trend.R) | confirmed |
| Fig 6A wider plotting helper | [pathogen_wider.R](../work_artifact/pathogen_wider.R) | helper |
| Fig 6A-B final aligned panels | [align_pathogen_pair.R](../work_artifact/align_pathogen_pair.R) | confirmed |
| Figure 1A: Baihe discharge source preparation | [prepare_and_verify.py](../outputs/baihe_monthly_indices_20260921/prepare_and_verify.py) | helper |
| Figure 1A: complete monthly grid preparation | [calculate_indices.R](../outputs/baihe_monthly_indices_20260921/calculate_indices.R) | helper |
| Figure 1A: calendar-month Gamma MLE SRI-1 and SDFAI | [calculate.R](../outputs/baihe_sri_matrix_20260922/calculate.R) | confirmed |
| Figure 1A: final SDFAI export | [export_plot.R](../outputs/baihe_gamma_sdfai_20260922/export_plot.R) | confirmed |
| Figure 1C-F: monthly HRT, TP, TN and chlorophyll a | [plot_monthly_environmental_curves_CNS4panel.R](../01_分析项目/环境变化曲线/plot_monthly_environmental_curves_CNS4panel.R) | confirmed |
| Figure 1B: monthly water level and net inflow | [plot_water_level_interannual_curve.R](../work/water_level_interannual_curve/plot_water_level_interannual_curve.R) | confirmed |
| Figure 1F source conversion | [export_source_csv.py](../work_artifact/chla_unit_fix_20260930/export_source_csv.py) | helper |
| Figure 1F vector export | [export_chla_svg.R](../work_artifact/chla_unit_fix_20260930/export_chla_svg.R) | helper |
| Figure 5B,C,E: water level, Chla and cyanobacteria CCF | [运行_CCF_两两关系.R](../01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/运行_CCF_两两关系.R) | confirmed |
| Figure 5D,F,G: TN:TP and temperature CCF | [补充_TOC_WT_原始水平CCF.R](../01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/补充_TOC_WT_原始水平CCF.R) | confirmed |
| Figure 5 original-level CCF panel plotting | [绘图_CCF_六关系_统一优化.R](../01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/绘图_CCF_六关系_统一优化.R) | confirmed |
| Figure 5 TN:TP and temperature panel plotting | [绘图_CCF_水温气温与整体蓝藻_统一单图.R](../01_分析项目/28_营养盐叶绿素蓝藻门CCF/分析代码/绘图_CCF_水温气温与整体蓝藻_统一单图.R) | confirmed |
| Cyanobacteria extraction configuration | [SEM_config.R](../01_分析项目/34_蓝藻门结构方程模型/配置/SEM_config.R) | helper |
| Figure 5A: extract sample-level cyanobacteria | [01_提取蓝藻门SEM每样本原始数据.R](../01_分析项目/34_蓝藻门结构方程模型/分析代码/01_提取蓝藻门SEM每样本原始数据.R) | helper |
| Figure 5A: monthly cyanobacterial abundance | [03_蓝藻门月均相对丰度纵向曲线.R](../01_分析项目/34_蓝藻门结构方程模型/分析代码/03_蓝藻门月均相对丰度纵向曲线.R) | confirmed |
| Figure 5A local launcher | [run_cyanobacteria_vertical_curve.R](../run_cyanobacteria_vertical_curve.R) | helper |
| Figure 5 final vector-panel layout | [layout_figure5_aligned.py](../work_artifact/layout_figure5_aligned.py) | helper |
| Figure 6B,C: ASV-level TITAN2 | [TITAN2_致病菌候选ASV_抽平数据.R](../01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_致病菌候选ASV_抽平数据.R) | confirmed |
| Figure 6B: ASV community sum(z) | [TITAN2_WL_ASV_群落sumz平滑曲线.R](../01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_WL_ASV_群落sumz平滑曲线.R) | confirmed |
| Figure 6C: top 20 unique-species representative ASVs | [TITAN2_WL_top20_不同种代表稳健ASV山脊图_最新注释.R](../01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_WL_top20_不同种代表稳健ASV山脊图_最新注释.R) | confirmed |
| Figure 6C: final unambiguous genus labels | [relabel_titan_genus_20260928.R](../work_artifact/relabel_titan_genus_20260928.R) | confirmed |
| ASV TITAN2 numerical QA | [TITAN2_WL_ASV_正式结果QA.R](../01_分析项目/36_TITAN2致病菌ASV水平_水位/分析代码/运行代码/TITAN2_WL_ASV_正式结果QA.R) | helper |
| Local TITAN ridge launcher | [replot_fig6b_titan_labels.R](../replot_fig6b_titan_labels.R) | helper |
| Supplementary Fig. 11: full richness and candidate fraction | [整体_Richness与候选致病菌相对丰度组合图.R](../01_分析项目/27_结果汇总/分析代码/运行代码/整体_Richness与候选致病菌相对丰度组合图.R) | confirmed |
| Supplementary Fig. 1: four-station study map | [plot_accurate_study_area.R](../figures/study_area/plot_accurate_study_area.R) | confirmed |
| GIS source inspection | [check_gis_source.R](../figures/study_area/check_gis_source.R) | helper |
