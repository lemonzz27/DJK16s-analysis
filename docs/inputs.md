# Input data

This is a code-only repository. Raw reads, full reference databases, source workbooks and cached model results are not uploaded with it. Use the same filenames and relative paths when staging the study data in a working copy. The local data catalogue is `provenance/input_files.tsv`; it lists source inputs and known working inputs, rather than guaranteeing that every dynamically constructed path has been captured.

Core inputs are the 36,097-ASV by 120-sample original count table, the whole-community table rarefied once to 21,000 reads per sample, representative sequences and rooted phylogeny, sample groups and station metadata, the environment table, SILVA taxonomy, the 2,755-candidate identifier set, and the MBPD filtered-hit and annotation tables. Candidate richness is calculated by subsetting the already rarefied complete-community table, without a second rarefaction. Candidate genus composition in the manuscript uses the original candidate abundance input.

The main input folders are `00_输入文件/01_ASV与序列`, `02_样本分组与站点`, `03_环境数据`, `04_分类注释`, `05_致病菌候选`, and `08_派生输入`. Names and columns remain those used by the original analyses. Original count matrices have ASVs in rows and sample IDs in columns. Group tables use sample/group fields; script-specific readers and the catalogue document other inputs.

Older analysis working folders under `99_运行环境/串行运行工作区/SILVA` contain renamed/staged inputs such as `tpm.csv`, `tax.csv`, `env.csv`, `groupsize.csv`, `groupmonth.csv` and `alpha1.csv`. They are data workspaces, not R package libraries. Several later corrections read these exact working copies, so merely supplying the core ASV table is insufficient for those scripts.

Hydrological preparation needs the historical Baihe monthly workbook and the 2024, 2025 and 2026 daily workbooks. The latter are mapped to `external_inputs/desktop/丹江口/磨峪湾白河流量数据(1)/成果2024-2026`. Study-area mapping needs the full shapefile sets, including .shx/.dbf/.prj/.cpg where present, with source and reuse terms recorded.

Final plotting scripts may consume cached GAMM predictions, TITAN tables, exported panels or other results produced earlier. Preserve these caches with their model-generation records. Their absence is an input dependency, not a reason to regenerate synthetic replacements. Manuscript figures assembled manually are not reconstructed automatically from a single command.

For a public release, deposit the permitted minimum dataset separately, add durable access identifiers and checksums, and supply restricted-data access instructions if necessary. Do not insert invented accession numbers or data DOIs.
