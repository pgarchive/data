# Code and data guide

## Primate Genome-phenome Archive: supplementary datasets and figure workflows

This guide accompanies the materials in `supplementary.information/`. It describes the supplementary datasets, the organization of the figure-analysis scripts and the current limitations affecting execution. Paths in this guide are relative to `supplementary.information/` unless stated otherwise. This document describes the files available at preparation time; it does not certify that all manuscript figures have been regenerated.

## 1. Overview

The supplementary materials comprise three datasets and the associated figure-analysis scripts. Data S1 contains phenotypic observations, trait metadata and group definitions. Data S2 contains genetic association and positional results. Data S3 contains functional enrichment and trait-pair overlap results. Each workbook includes an index of contents.

Main-figure scripts are located in `workflow/MAIN/Scripts/`. Individual supplementary-figure scripts are located in `workflow/SUPP/S1/` through `S50/`; the S11 directory is named `s11` with a lowercase s. Six additional scripts are stored in `workflow/SUPP/S6,9,12,15,18,21_MERGED/`. Script filenames identify the associated figures or panels, while their contents document data loading, filtering, joins, statistical calculations and plotting choices.

The supplied collection retains some paths from the original analysis environment and seven scripts with syntax errors. These limitations are documented in Sections 4 and 5. The final inventory records the relative path, MD5 checksum and a brief description of every file under `tables/` and `workflow/`.

## 2. Supplementary datasets

The `tables/` directory contains three original Excel workbooks (`.xlsx`) and three OpenDocument Spreadsheet exports (`.ods`) with matching basenames. The ODS files provide open-format copies of the same datasets. They are alternatives to the Excel files, not additional observations. Open one format per dataset to avoid double counting.

### Data S1: phenotypic data and trait metadata

`tables/Data_S1.xlsx` and `tables/Data_S1.ods` contain nine worksheets: an index of contents; non-human primate phenomic data; trait classification and data sources; human-equivalent measurements; quantitative and qualitative trait metadata; trait redundancy; an allometric trait list; and trait top/bottom selections. Original worksheet names, including abbreviations and truncations, are preserved.

### Data S2: genetic results and positional annotations

`tables/Data_S2.xlsx` and `tables/Data_S2.ods` contain seven worksheets: an index of contents; filtered gene-alignment metadata; CAAS contrast metadata; CAAS results for traits with contrasts; an RER results table labeled FDR < 0.05; CAAS mutation results for traits with contrasts; and mammalian-equivalent positions. CAAS denotes convergent amino acid substitutions; RER denotes relative evolutionary rate.

### Data S3: enrichment and overlap results

`tables/Data_S3.xlsx` and `tables/Data_S3.ods` contain five worksheets: an index of contents; two Permul functional-enrichment tables; a PGLS functional-enrichment table; and trait-pair overlap results. Worksheet labels identify the FDR filtering used where specified. PGLS denotes phylogenetic generalized least squares.

### Conversion checks

The ODS exports were checked against all 21 original worksheets, covering 1,481,536 non-empty cells. Cell contents and worksheet names/order were compared, accounting for numeric precision, ODS whitespace encoding and an Excel date serial. No substantive cell-content differences were detected. The source workbooks contain no cell formulas in the inspected worksheet XML. This comparison verifies tabular contents; it is not a claim of identical rendering across spreadsheet applications or validation of the underlying scientific analyses.

## 3. Workflow organization and file formats

There are 69 R scripts: 13 in `workflow/MAIN/Scripts/`, 50 individual supplementary-figure scripts, and six merged scripts. The supplementary folders also contain auxiliary datasets, configuration files and serialized R objects.

The main shared input directory is `workflow/MAIN/InputFiles/`. Inputs are also stored in individual supplementary-figure folders, including `SUPP/S27/Configs_raw/`, `SUPP/S36/PHYLOP_files/` and `SUPP/S37/VARIATION_files/`. There is no single shared `input/` directory in the supplied structure. Scripts may depend on data supplied under another figure's folder.

Plain-text inputs use extensions such as `.txt`, `.tsv`, `.tab` and `.out`; some text tables have no extension. Inspect the script's separator and header arguments before importing a table: extensions alone do not establish whether a header is present. Newick files contain phylogenetic trees. Files read with `readRDS()` are serialized R objects, including plotting objects and color palettes. Several have a `.gz` suffix. Do not interpret every `.gz` file as a compressed text table; follow the reader used in the corresponding script.

For example, the available trait palette can be inspected in R from the repository root with `readRDS("supplementary.information/workflow/MAIN/InputFiles/traits_palette.gz")`. Loading an object does not run the associated analysis or produce a figure.

Many scripts save R objects through `saveRDS()` rather than exporting a final image. Such objects are outputs, even where the filename has no extension. Check the actual saving commands when evaluating what a script produces. Calls to `View()` are intended for interactive inspection and may require adjustment for non-interactive execution.

## 4. Input resolution and execution requirements

Relative paths in R are resolved against the working directory, not automatically against the script's location. Running `Rscript` with an absolute script path does not itself change this behavior. Work on a copy if adapting paths for review, and record any changes and the R/package versions used.

### Known path mismatches

- Many supplementary scripts reference `../../Figures/Fig1/...` or `../../Figures/Fig2/...`. These directories are absent from the supplied workflow layout. Candidate inputs often exist in `workflow/MAIN/InputFiles/` or another supplementary folder.
- Main scripts frequently use bare filenames while residing in `MAIN/Scripts/`, separate from `MAIN/InputFiles/`. Choosing `MAIN/InputFiles/` as the working directory may resolve some references, but does not fix every legacy path such as `../Fig1/...`.
- Several `readRDS()` calls omit the `.gz` suffix of the available object, for example `traits_palette`, `joined_LOWESTPVAL_FEB24` and the `PHYLOP_files/SCENARIO1` object. R requires the actual filename; compression support does not make suffixes optional.
- `MAIN/Scripts/Figure2F_2_COUNTS.R` requests `traits_palette.rds`, whereas the available candidate is `MAIN/InputFiles/traits_palette.gz`. Confirm the expected object before substituting it.
- Repeated basenames can occur in different folders. A matching basename is a candidate location, not proof that two files are interchangeable. Use the per-file checksums to distinguish copies.

A static scan found at least one unresolved input reference in 68 of the 69 scripts when relative paths were interpreted from each script's directory. This is a diagnostic of the current layout, not a result from executing the analyses. Dynamic inputs, intermediate objects and branch-specific code can require additional checks. No universal working directory or launch command has been validated for the entire collection.

### Referenced inputs not located in the supplied workflow

- `SUPP/S15/FigS15.R`: `Metadata_numeric_traits_2.0.tsv`, `Metadata_categorical_traits_FINAL.tsv` and `Metadata_binary_traits_FINAL.tsv` are referenced through legacy results paths. A separate absolute path requests `FINAL_TREE_BL.newick`. A differently named tree is supplied, but equivalence has not been established.
- `SUPP/S26/FigS26.R`: `all_categorical_traits_extremes_2BG_FINAL` is requested through `readRDS()` and was not found.
- `SUPP/S37/FigS37.R`: `ALL_CAAS_PHYLOPs_APR23.txt` is requested under `S36/PHYLOP_files/caas_phyloP/` and was not found. Chromosome-specific files exist there, but substituting or combining them has not been validated.

These observations concern the supplied `tables/` and `workflow/` materials. Likewise, `Cluster_biometabolites` is read by S15, but is written earlier in that script and should be treated as an intermediate output rather than automatically as a missing source file.

### Software environment

Execution requires R and the packages loaded by each selected script. Examples include ggplot2, dplyr, tidyr, ape, phytools, ComplexHeatmap, FactoMineR, ggtree and WebGestaltR. This is not a complete installation list. Some scripts include package-installation or external-service calls; inspect these before running them. The supplied workflow does not include a dependency lockfile or a validated environment specification. Package versions and full end-to-end figure regeneration have not been validated in this review.

## 5. Warning: retained syntax errors in seven scripts

**Seven scripts contain paired triple double-quote delimiters (`"""`). These delimiters are not valid R block comments. They cause parsing errors and prevent the complete file from being executed with `Rscript` or `source()`. The problem does not merely deactivate the enclosed code.**

The affected files and delimiter lines are:

- `workflow/MAIN/Scripts/Fig5C.R`: lines 217 and 373, around additional selection/CAAS overlap analyses.
- `workflow/SUPP/S15/FigS15.R`: lines 823 and 834, around assignments of redundant-trait groups.
- In `workflow/SUPP/S6,9,12,15,18,21_MERGED/`, `S6_MERGED.R`, `S12_MERGED.R`, `S15_MERGED.R` and `S18_MERGED.R`: lines 263 and 274, around a `TestCor` / `ApplyFdrCor()` block.
- In the same directory, `S9_MERGED.R`: lines 262 and 273, around the corresponding `TestCor` block. `S21_MERGED.R` is not among the seven files with this error.

These delimiters have been retained because it is not established whether the enclosed code was intended to run. Removing only the delimiters would activate the enclosed analyses; commenting every line of the block would leave it inactive. Neither interpretation has been adopted in the supplied seven files. Readers can inspect the code, but execution requires clarification and an explicit correction.

A separate trailing `+` after `saveRDS(p2, "FigureS26")` was removed from `workflow/SUPP/S26/FigS26.R`. That script now passes R's syntax parser. A fresh syntax check of all 69 files found 62 that parse and the seven failures listed above. Passing the parser does not establish that inputs, packages, object definitions or scientific results are correct. No other script corrections are represented by this guide.

## 6. Integrity checks and scope of the inventory

The inventory below lists every regular file recursively under `tables/` and `workflow/`, including auxiliary inputs and all original/exported workbook pairs. It contains 356 entries: six table files and 350 workflow files. This guide and its PDF export are outside those two directories and are therefore excluded, avoiding a self-referential checksum.

MD5 values were computed from each file's complete stored bytes, including the compression wrapper for `.gz` files. They identify this snapshot and can detect changed or incomplete downloads; they are not checksums of decompressed content. Editing a file, re-saving a workbook or recompressing an object can change its MD5 without necessarily changing the scientific data. Descriptions summarize the visible filename, file content or its use in the scripts; they do not substitute for a formal data dictionary or establish undocumented scenario meanings.

To check an individual file, run the following from `supplementary.information/` and compare the returned digest with its inventory entry:

```bash
md5sum "tables/Data_S1.ods"
```

On macOS, the corresponding command is `md5 "tables/Data_S1.ods"`.

## 7. Complete file inventory

Entries are grouped by directory and sorted by relative path. Each entry gives the exact filename, its MD5 digest and a short description.

### tables

- **File:** `tables/Data_S1.ods`  
  **MD5:** `c7654982c600f91f79364c66ce571017`  
  **Description:** Data S1: phenotypic observations, trait metadata and extreme selections; OpenDocument Spreadsheet export.

- **File:** `tables/Data_S1.xlsx`  
  **MD5:** `56ea80211386017989f689bc45bcf309`  
  **Description:** Data S1: phenotypic observations, trait metadata and extreme selections; original Excel workbook.

- **File:** `tables/Data_S2.ods`  
  **MD5:** `e49c8bd2942568d150e70ed1ccc316c6`  
  **Description:** Data S2: genetic results, alignment metadata and positional annotations; OpenDocument Spreadsheet export.

- **File:** `tables/Data_S2.xlsx`  
  **MD5:** `c47e3d982f065f4bca50d5a1866c2270`  
  **Description:** Data S2: genetic results, alignment metadata and positional annotations; original Excel workbook.

- **File:** `tables/Data_S3.ods`  
  **MD5:** `2be36814b2e568ec971758d70e9593b1`  
  **Description:** Data S3: functional enrichment and trait-pair overlaps; OpenDocument Spreadsheet export.

- **File:** `tables/Data_S3.xlsx`  
  **MD5:** `eab1fbe0e04e6c6ed21724e0e52e1948`  
  **Description:** Data S3: functional enrichment and trait-pair overlaps; original Excel workbook.

### workflow/MAIN/InputFiles

- **File:** `workflow/MAIN/InputFiles/ALL_BM_OCT23_BETAS.txt`  
  **MD5:** `5c26f5b21299718412fad57a63d784de`  
  **Description:** Brownian-motion model association results, including coefficient information.

- **File:** `workflow/MAIN/InputFiles/ALL_ENRICHMENTS_SEP23_FINAL.txt`  
  **MD5:** `34da72cde81b459bf4e55afd793f3ff7`  
  **Description:** Functional enrichment results used in figure scripts.

- **File:** `workflow/MAIN/InputFiles/ALL_LAMBDA_OCT23_BETAS.txt`  
  **MD5:** `668ccbd0916505ca29551e41999e98ff`  
  **Description:** Lambda model association results, including coefficient information.

- **File:** `workflow/MAIN/InputFiles/ALL_RER_CONVERGE_UNIQUE_MAR23.txt`  
  **MD5:** `42b6953d8f7d4f6754239ebad375b2ef`  
  **Description:** RERconverge results used by association figure scripts.

- **File:** `workflow/MAIN/InputFiles/ALL_SELECTION_STUDY_CASES_FILT.txt`  
  **MD5:** `2ad090469a2ff732d90c9dfd4a296531`  
  **Description:** Filtered selection results for the illustrated study cases.

- **File:** `workflow/MAIN/InputFiles/CAAS_FINAL_MAR23.tab`  
  **MD5:** `b17b90c51fc92da5da8a222748f04d8f`  
  **Description:** CAAS results used by association and comparative figure scripts.

- **File:** `workflow/MAIN/InputFiles/CAAS_global_and_altern.txt`  
  **MD5:** `f5c100a0c1193c866e89ea286404503c`  
  **Description:** CAAS results comparing global and alternative group configurations.

- **File:** `workflow/MAIN/InputFiles/DATASET_CAAS_OVERLAPS_MAR24.gz`  
  **MD5:** `9925b48171d0b01b464d3779787b5f4a`  
  **Description:** Compressed serialized R dataset of CAAS overlaps.

- **File:** `workflow/MAIN/InputFiles/DATASET_RERS_OVERLAPS_MAR24.gz`  
  **MD5:** `6172d2817c6d323c5db0f011a580084d`  
  **Description:** Compressed serialized R dataset of RER overlaps.

- **File:** `workflow/MAIN/InputFiles/FINAL_TREE_BL_WITHOUTNODELAB.newick`  
  **MD5:** `d9a6cb33c19f84a5335a42f958558735`  
  **Description:** Phylogenetic tree in Newick format used by primate figure scripts.

- **File:** `workflow/MAIN/InputFiles/Maximum.lifespan..y._altern1_config.txt`  
  **MD5:** `456c08980f51186f4d243a6617216341`  
  **Description:** Species/group configuration for Maximum lifespan y altern1; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/Maximum.lifespan..y._altern2_config.txt`  
  **MD5:** `032b91cfa4ffc17b4ae9dbcdd658c759`  
  **Description:** Species/group configuration for Maximum lifespan y altern2; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/Maximum_lifespan__y__config.txt`  
  **MD5:** `96e917c2f40710dc9ea7abd39343916f`  
  **Description:** Species/group configuration for Maximum lifespan y; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/PRIMATES_PHENOMIC_DATASET.tsv`  
  **MD5:** `2547224c5fc0845fa6e7a17bd6cf166e`  
  **Description:** Non-human primate phenotypic observations used by the figure scripts.

- **File:** `workflow/MAIN/InputFiles/Qualitative_extremes`  
  **MD5:** `b35618c7f0b87cb77d7c69db07bcd7f7`  
  **Description:** Qualitative trait group and extreme-selection table.

- **File:** `workflow/MAIN/InputFiles/Quantitative_extremes`  
  **MD5:** `c51a85a75150075e2e231e643ba4219a`  
  **Description:** Quantitative trait extreme-selection table.

- **File:** `workflow/MAIN/InputFiles/Traits_clusters_FINAL_JAN23.tsv`  
  **MD5:** `b22881816aa5f164be99769ac02a42ae`  
  **Description:** Trait clustering assignments used in coverage and clustering figures.

- **File:** `workflow/MAIN/InputFiles/Traits_domains_FIG2.txt`  
  **MD5:** `60d09ffdbf0c229e6d92d4c02184f0f1`  
  **Description:** Trait domain assignments used by Figure 2 scripts.

- **File:** `workflow/MAIN/InputFiles/Traits_domains_SPLIT_FEB22.txt`  
  **MD5:** `60d09ffdbf0c229e6d92d4c02184f0f1`  
  **Description:** Trait domain and secondary-domain annotations.

- **File:** `workflow/MAIN/InputFiles/Traits_domains_SPLIT_FEB22_ALL_UNDERSCORES.txt`  
  **MD5:** `b13cd8aabed6f3440100f6c6bd04b347`  
  **Description:** Trait-domain annotations with underscore-based identifiers.

- **File:** `workflow/MAIN/InputFiles/Traits_domains_sorted_FIG1.txt`  
  **MD5:** `6aa6c096d1dced5236279daac75bce7c`  
  **Description:** Ordered trait domain annotations used by Figure 1 and related scripts.

- **File:** `workflow/MAIN/InputFiles/TrophicGuild_Insectivore_altern1_config.txt`  
  **MD5:** `37e49bdf95afb2db53fb2696c8c4a25c`  
  **Description:** Species/group configuration for TrophicGuild Insectivore altern1; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/TrophicGuild_Insectivore_altern2_config.txt`  
  **MD5:** `4ba93385ef89742fb6f9f39eef2023ee`  
  **Description:** Species/group configuration for TrophicGuild Insectivore altern2; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/TrophicGuild_Insectivore_config.txt`  
  **MD5:** `cd188f92915a9d680559a3ee7c84285a`  
  **Description:** Species/group configuration for TrophicGuild Insectivore; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/White.Blood.Cells_altern1_config.txt`  
  **MD5:** `fb5f7abd06d5864a7c3aa7c39c807f6a`  
  **Description:** Species/group configuration for White Blood Cells altern1; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/White.Blood.Cells_altern2_config.txt`  
  **MD5:** `60c042ddf40797ab0ca58fdb5cbb987e`  
  **Description:** Species/group configuration for White Blood Cells altern2; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/White.Blood.Cells_config.txt`  
  **MD5:** `2aec7f968dd5476620ea72d3ff714bf1`  
  **Description:** Species/group configuration for White Blood Cells; used to define or compare trait groups.

- **File:** `workflow/MAIN/InputFiles/genes_coordinates_sorted.txt`  
  **MD5:** `ef4e360f8196f579ff60159267255b3a`  
  **Description:** Sorted gene coordinates used in Manhattan plots.

- **File:** `workflow/MAIN/InputFiles/gg3_FEB24_contr.gz`  
  **MD5:** `bd23e2ebaab0cfebe27f080f80021da5`  
  **Description:** Compressed serialized R object used in contrast-overlap plotting.

- **File:** `workflow/MAIN/InputFiles/human_direct_measurements.tsv`  
  **MD5:** `96f72fdd09dda96487b3ca6fe8b72a9c`  
  **Description:** Human measurements used for comparison with primate traits.

- **File:** `workflow/MAIN/InputFiles/joined_LOWESTPVAL_FEB24.gz`  
  **MD5:** `d41e39cdb8146496842166e7c6514ddf`  
  **Description:** Compressed serialized R association data used in Manhattan plots.

- **File:** `workflow/MAIN/InputFiles/joined_overlap_FEB24.txt`  
  **MD5:** `1c08d6c542eb5a0388c03f4f35611263`  
  **Description:** Joined overlap results used by figure scripts.

- **File:** `workflow/MAIN/InputFiles/names_FINAL.txt`  
  **MD5:** `736179c8c9a8007cc407095a523c8b97`  
  **Description:** Species-name list used to identify sequenced primates.

- **File:** `workflow/MAIN/InputFiles/primate_families_palette.rds`  
  **MD5:** `228a4769000d2e195d133cdc3d8e8114`  
  **Description:** Serialized R color palette for primate families.

- **File:** `workflow/MAIN/InputFiles/traits_palette.gz`  
  **MD5:** `3e869d6cf6727abcd099d45dc59848b2`  
  **Description:** Compressed serialized R color palette for trait groups.

### workflow/MAIN/Scripts

- **File:** `workflow/MAIN/Scripts/Fig5A2.R`  
  **MD5:** `704eacdbfb926d3544c65d21a4df6294`  
  **Description:** Maximum-lifespan global/alternative group comparison for Figure 5A.

- **File:** `workflow/MAIN/Scripts/Fig5A3.R`  
  **MD5:** `397a811d7ff7d4fad7496fe2860f0039`  
  **Description:** White-blood-cell global/alternative group comparison for Figure 5A.

- **File:** `workflow/MAIN/Scripts/Fig5B.R`  
  **MD5:** `456ff34cbf59be48fc50c16cf89ae467`  
  **Description:** Global/alternative CAAS group comparisons for Figure 5B.

- **File:** `workflow/MAIN/Scripts/Fig5C.R`  
  **MD5:** `6892375fb3c5ac969258fc00470c48b2`  
  **Description:** Selection and CAAS overlap analysis for Figure 5C; retained syntax error.

- **File:** `workflow/MAIN/Scripts/Figure1A.R`  
  **MD5:** `f0279314f11bb96c60be1d89e09c8677`  
  **Description:** Trait coverage heatmap and clustering visualization for Figure 1A.

- **File:** `workflow/MAIN/Scripts/Figure1B.R`  
  **MD5:** `2c236f7230bcba40d97073cd4ff4be9c`  
  **Description:** Extreme-group contrast summaries by primate family for Figure 1B.

- **File:** `workflow/MAIN/Scripts/Figure1C1_STRIPCHART.R`  
  **MD5:** `dce5df7b8ae73b40d0c145226ce8ead6`  
  **Description:** Body-mass extreme-selection strip chart for Figure 1C.

- **File:** `workflow/MAIN/Scripts/Figure1C2_STRIPCHART.R`  
  **MD5:** `95268629d7d7860c2c87f13459f867e3`  
  **Description:** Lifespan and qualitative extreme-selection plots for Figure 1C.

- **File:** `workflow/MAIN/Scripts/Figure2A.R`  
  **MD5:** `56e61c878c0d209d6584c3270865d93c`  
  **Description:** FDR-based association counts across CAAS, RER and model results for Figure 2A.

- **File:** `workflow/MAIN/Scripts/Figure2B_2_LOWEST.R`  
  **MD5:** `70947d8d76f3819ec0fd35e694c41620`  
  **Description:** Manhattan plots of association results for Figure 2B.

- **File:** `workflow/MAIN/Scripts/Figure2D_E.R`  
  **MD5:** `bb596656ffc9c9529008a6732fb0e731`  
  **Description:** Functional enrichment visualization for Figure 2D/E.

- **File:** `workflow/MAIN/Scripts/Figure2F_2_COUNTS.R`  
  **MD5:** `8daa6c48675ee6ebfe393bd5db16540d`  
  **Description:** Overlap matrices and counts for Figure 2F.

- **File:** `workflow/MAIN/Scripts/Figure5A.R`  
  **MD5:** `b963bac66443826617aea132bdbe7b50`  
  **Description:** Phylogenetic and trait-group visualization for Figure 5A.

### workflow/SUPP/S1

- **File:** `workflow/SUPP/S1/FigS1.R`  
  **MD5:** `7b5542ef51981473271670db41e31f47`  
  **Description:** R analysis and plotting script for Supplementary Figure 1.

### workflow/SUPP/S10

- **File:** `workflow/SUPP/S10/FigS10.R`  
  **MD5:** `8881dfebba832f8ad7473a14232a9083`  
  **Description:** R analysis and plotting script for Supplementary Figure 10.

### workflow/SUPP/S12

- **File:** `workflow/SUPP/S12/FigS12.R`  
  **MD5:** `5fc26cd923878b3770e01267376c878d`  
  **Description:** R analysis and plotting script for Supplementary Figure 12.

### workflow/SUPP/S13

- **File:** `workflow/SUPP/S13/FigS13.R`  
  **MD5:** `5c9314dfc168289c7096e794b69bfd17`  
  **Description:** R analysis and plotting script for Supplementary Figure 13.

### workflow/SUPP/S14

- **File:** `workflow/SUPP/S14/FigS14.R`  
  **MD5:** `13c01d540dbd10f09c7e97def7975d2a`  
  **Description:** R analysis and plotting script for Supplementary Figure 14.

### workflow/SUPP/S15

- **File:** `workflow/SUPP/S15/FigS15.R`  
  **MD5:** `895d4919a88b319daaafdd1fad77b6c7`  
  **Description:** R analysis and plotting script for Supplementary Figure 15; retained syntax error.

### workflow/SUPP/S16

- **File:** `workflow/SUPP/S16/FigS16.R`  
  **MD5:** `94bd79178fc976efd987afed19d4cc75`  
  **Description:** R analysis and plotting script for Supplementary Figure 16.

### workflow/SUPP/S17

- **File:** `workflow/SUPP/S17/FigS17.R`  
  **MD5:** `3e65e87aeec04da8e0de0b3b5155fee4`  
  **Description:** R analysis and plotting script for Supplementary Figure 17.

### workflow/SUPP/S18

- **File:** `workflow/SUPP/S18/FigS18.R`  
  **MD5:** `b1cdf118f045238827199b62f5d8f13c`  
  **Description:** R analysis and plotting script for Supplementary Figure 18.

### workflow/SUPP/S19

- **File:** `workflow/SUPP/S19/FigS19.R`  
  **MD5:** `d8d2b3d12df7061b0b064a3ba21b3470`  
  **Description:** R analysis and plotting script for Supplementary Figure 19.

### workflow/SUPP/S2

- **File:** `workflow/SUPP/S2/FigS2.R`  
  **MD5:** `c882ce95e91f5516b780e14aaa748bbe`  
  **Description:** R analysis and plotting script for Supplementary Figure 2.

### workflow/SUPP/S20

- **File:** `workflow/SUPP/S20/FigS20.R`  
  **MD5:** `871f04044bdde3f9948877a5ee1455af`  
  **Description:** R analysis and plotting script for Supplementary Figure 20.

- **File:** `workflow/SUPP/S20/all_cors_pasted_JAN23.txt`  
  **MD5:** `5e90b52dbfc8f50c81c533ba3618c015`  
  **Description:** Combined trait-correlation results for Supplementary Figure 20.

- **File:** `workflow/SUPP/S20/qual_dummy_JAN23`  
  **MD5:** `a2fa28267aa8b9cd985fd9afbf89ea28`  
  **Description:** Dummy-coded qualitative trait data for Supplementary Figure 20.

### workflow/SUPP/S21

- **File:** `workflow/SUPP/S21/FigS21.R`  
  **MD5:** `711269be3e7ba0850b2cdc3348fe4533`  
  **Description:** R analysis and plotting script for Supplementary Figure 21.

- **File:** `workflow/SUPP/S21/pruned_tree.newick`  
  **MD5:** `f476d99169cbb4229b69b42c1e9cfcb4`  
  **Description:** Phylogenetic tree in Newick format pruned for Supplementary Figure 21.

### workflow/SUPP/S22

- **File:** `workflow/SUPP/S22/AAs_colors.gz`  
  **MD5:** `141f8c5da8d5c18d697131abda5d722f`  
  **Description:** Compressed serialized R amino-acid color palette.

- **File:** `workflow/SUPP/S22/Allometric_traits.txt`  
  **MD5:** `2f251392842d370c218347cdc4f79cd0`  
  **Description:** Allometric trait data used by Supplementary Figure 22.

- **File:** `workflow/SUPP/S22/FigS22.R`  
  **MD5:** `38b95f9a1136f69d6186e37e8409002b`  
  **Description:** R analysis and plotting script for Supplementary Figure 22.

- **File:** `workflow/SUPP/S22/Metadata_continous_traits_FINAL.tsv`  
  **MD5:** `e2dbde3de3214093357a45709b11b049`  
  **Description:** Continuous-trait metadata used by Supplementary Figure 22.

### workflow/SUPP/S23

- **File:** `workflow/SUPP/S23/FigS23.R`  
  **MD5:** `4fbef666c4f2770d52c06df635cac97a`  
  **Description:** R analysis and plotting script for Supplementary Figure 23.

- **File:** `workflow/SUPP/S23/n_species_df.gz`  
  **MD5:** `e682c32047d54e433efc4357f3d3e0b6`  
  **Description:** Compressed serialized R species-count object for Supplementary Figure 23.

### workflow/SUPP/S24

- **File:** `workflow/SUPP/S24/FigS24.R`  
  **MD5:** `2337cd8f4bfd3459bfe9e9b9de543f87`  
  **Description:** R analysis and plotting script for Supplementary Figure 24.

- **File:** `workflow/SUPP/S24/NONALLOMETRIC_EXTREMES`  
  **MD5:** `ef7a0b96babef691499af7a0dabcf082`  
  **Description:** Non-allometric trait extreme selections for Supplementary Figure 24.

### workflow/SUPP/S25

- **File:** `workflow/SUPP/S25/ALLOMETRIC_EXTREMES`  
  **MD5:** `c44976e4c5206e78e603c705fe0371b5`  
  **Description:** Allometric trait extreme selections for Supplementary Figure 25.

- **File:** `workflow/SUPP/S25/FigS25.R`  
  **MD5:** `93b5cc58b747dc6994efe489741e9893`  
  **Description:** R analysis and plotting script for Supplementary Figure 25.

### workflow/SUPP/S26

- **File:** `workflow/SUPP/S26/FigS26.R`  
  **MD5:** `1360fad87bc2911872585ece8ee00341`  
  **Description:** R analysis and plotting script for Supplementary Figure 26; trailing syntax error corrected.

### workflow/SUPP/S27/Configs_raw

- **File:** `workflow/SUPP/S27/Configs_raw/Basophils_Abs_config.txt`  
  **MD5:** `60ca8cd6ddf551e42001e2c6ad31277f`  
  **Description:** Species/group configuration for Basophils Abs; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Basophils___config.txt`  
  **MD5:** `928e165e5786e41112e4a1801a5f436e`  
  **Description:** Species/group configuration for Basophils; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Blood_Urea_Nitrogen_config.txt`  
  **MD5:** `9a286e399dea387328b5b1e1139cae16`  
  **Description:** Species/group configuration for Blood Urea Nitrogen; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Bun_Creat_Ratio_config.txt`  
  **MD5:** `8ca5c79ae2657bdf211746f90981a0bf`  
  **Description:** Species/group configuration for Bun Creat Ratio; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Cholesterol_config.txt`  
  **MD5:** `fd1d205c6c48bc63b3a200e65e2be82f`  
  **Description:** Species/group configuration for Cholesterol; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Diencephalon_config.txt`  
  **MD5:** `5125952b47bb03a8655b0fd627760261`  
  **Description:** Species/group configuration for Diencephalon; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/GGT_config.txt`  
  **MD5:** `ac88c3b61796f56c0d82cbf0f7157238`  
  **Description:** Species/group configuration for GGT; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Globulin_config.txt`  
  **MD5:** `ee6b789fea2c6923a2c039ca92c7f312`  
  **Description:** Species/group configuration for Globulin; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Hippocampus_config.txt`  
  **MD5:** `431ed57c33915f08b6af24d77b687f0f`  
  **Description:** Species/group configuration for Hippocampus; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Hypothalamus_config.txt`  
  **MD5:** `d5bad1920a97b25f739f57b13ca8ed8e`  
  **Description:** Species/group configuration for Hypothalamus; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/MOB_config.txt`  
  **MD5:** `69ac2bf87662ceb143111b0633fddddd`  
  **Description:** Species/group configuration for MOB; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/M_Max_M_NeonateWeight_g_config.txt`  
  **MD5:** `5264e3e79dd850a295a76c01b607e52d`  
  **Description:** Species/group configuration for M Max M NeonateWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Maximum_lifespan__y__config.txt`  
  **MD5:** `96e917c2f40710dc9ea7abd39343916f`  
  **Description:** Species/group configuration for Maximum lifespan y; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_AdiposeDepots_g_config.txt`  
  **MD5:** `356430e1b33a4d4232c181eb254d0b3c`  
  **Description:** Species/group configuration for Phylo Permul AdiposeDepots g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_AdultFemale_config.txt`  
  **MD5:** `398beabb97aba24aeb2f79fe624c8553`  
  **Description:** Species/group configuration for Phylo Permul AdultFemale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_AdultMales_config.txt`  
  **MD5:** `531773aaf5a1d80c6b6bcd7b05da3acb`  
  **Description:** Species/group configuration for Phylo Permul AdultMales; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Age_of_female_at_first_reproduction__y__config.txt`  
  **MD5:** `0ed5867880bba9eda34b2dc640829d27`  
  **Description:** Species/group configuration for Phylo Permul Age of female at first reproduction y; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_AgilityScore_config.txt`  
  **MD5:** `f347574834fb6b1c4165c3863f40f51c`  
  **Description:** Species/group configuration for Phylo Permul AgilityScore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Alkaline_Phosphatase_config.txt`  
  **MD5:** `75ed1e0f08f1bdc26ab788a6431bbbef`  
  **Description:** Species/group configuration for Phylo Permul Alkaline Phosphatase; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Amygdala_config.txt`  
  **MD5:** `ef3dc9c99a75b73a114eb9703fcf71e1`  
  **Description:** Species/group configuration for Phylo Permul Amygdala; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Bilirubin_Total_config.txt`  
  **MD5:** `eaf796c5a1e91b5816663eebd9f4e1b3`  
  **Description:** Species/group configuration for Phylo Permul Bilirubin Total; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_BodyHair_config.txt`  
  **MD5:** `432af3afde84f3035ce1f019ae67a0b5`  
  **Description:** Species/group configuration for Phylo Permul BodyHair; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_BodyMassFemale_kg_config.txt`  
  **MD5:** `266a725282ffb1195d242212eb0395f7`  
  **Description:** Species/group configuration for Phylo Permul BodyMassFemale kg; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_BodyMassMale_kg_config.txt`  
  **MD5:** `a01384dc2bb5f9c53557b3d670fc3ae3`  
  **Description:** Species/group configuration for Phylo Permul BodyMassMale kg; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_BodyMass_kg_config.txt`  
  **MD5:** `3627bd9d3d99c626ec706a5655c009a9`  
  **Description:** Species/group configuration for Phylo Permul BodyMass kg; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_BodyPelageDichrom_config.txt`  
  **MD5:** `378cd280030eca5b950c2f6719c23e3f`  
  **Description:** Species/group configuration for Phylo Permul BodyPelageDichrom; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_CPK_config.txt`  
  **MD5:** `ed79e07ccadd2e8e223b99a4d5f74212`  
  **Description:** Species/group configuration for Phylo Permul CPK; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Coenspecific_violence_config.txt`  
  **MD5:** `04fb8a706f8fdcc54e4a464984393d23`  
  **Description:** Species/group configuration for Phylo Permul Coenspecific violence; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Colour_vision_Dichromatic_config.txt`  
  **MD5:** `899f495d7120c877343f370e006562dd`  
  **Description:** Species/group configuration for Phylo Permul Colour vision Dichromatic; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Colour_vision_Monochromatic_config.txt`  
  **MD5:** `9af1878c4cae5eb5c19de16e7fe6bb4a`  
  **Description:** Species/group configuration for Phylo Permul Colour vision Monochromatic; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Colour_vision_Polymorphic_config.txt`  
  **MD5:** `0da27ffa1627cdfe7f0433b24577c3ce`  
  **Description:** Species/group configuration for Phylo Permul Colour vision Polymorphic; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Colour_vision_Trichromatic_config.txt`  
  **MD5:** `41680de43ae801470ca63f4a06baca21`  
  **Description:** Species/group configuration for Phylo Permul Colour vision Trichromatic; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_CombinedTestesMass_in_g_config.txt`  
  **MD5:** `d3a94222fbac6c5c40001b8d0a56b2a8`  
  **Description:** Species/group configuration for Phylo Permul CombinedTestesMass in g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_DayLength_km_config.txt`  
  **MD5:** `a4ed7ddbeaff1d06761b3dfe26981d1e`  
  **Description:** Species/group configuration for Phylo Permul DayLength km; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_DielActivity_CA_config.txt`  
  **MD5:** `5b5a9b58594e9297fd1e68e5ef89e9ff`  
  **Description:** Species/group configuration for Phylo Permul DielActivity CA; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_DielActivity_D_config.txt`  
  **MD5:** `1c35c958e65d1589c47ab09194d7a0ed`  
  **Description:** Species/group configuration for Phylo Permul DielActivity D; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_DielActivity_N_config.txt`  
  **MD5:** `143cfec045201f4bee60acdd6f950daa`  
  **Description:** Species/group configuration for Phylo Permul DielActivity N; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Dorsal_Luminance_config.txt`  
  **MD5:** `d867e13dc19461b52a6c37a3813180fb`  
  **Description:** Species/group configuration for Phylo Permul Dorsal Luminance; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Eosinophils_Abs_config.txt`  
  **MD5:** `6ae11fc2b62c86104f1b0db42b65781e`  
  **Description:** Species/group configuration for Phylo Permul Eosinophils Abs; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Eosinophils___config.txt`  
  **MD5:** `53eb61b4eee45e211ad6100732af2c08`  
  **Description:** Species/group configuration for Phylo Permul Eosinophils; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Fauna_config.txt`  
  **MD5:** `087a0078f27b5e2d69349730f0c3303d`  
  **Description:** Species/group configuration for Phylo Permul Fauna; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Female_Ranging_patterns_config.txt`  
  **MD5:** `2dd4a61fd4b03a1a1332d94d9dedc6e3`  
  **Description:** Species/group configuration for Phylo Permul Female Ranging patterns; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_GR_MidRangeLat_dd_config.txt`  
  **MD5:** `f2c4afd9d00678e32a1230a195b7f26b`  
  **Description:** Species/group configuration for Phylo Permul GR MidRangeLat dd; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_GenitalSkin_config.txt`  
  **MD5:** `183d0ce768bd0d67d73187ac08faa4ee`  
  **Description:** Species/group configuration for Phylo Permul GenitalSkin; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Gestation__days__config.txt`  
  **MD5:** `e7f3e32e29fc77761ab369a97f02f5bd`  
  **Description:** Species/group configuration for Phylo Permul Gestation days; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Habitat_Forest_Savanna_Shrubland_config.txt`  
  **MD5:** `eb4da763485021a8e792af736a7d7422`  
  **Description:** Species/group configuration for Phylo Permul Habitat Forest Savanna Shrubland; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Habitat_Forest_Savanna_config.txt`  
  **MD5:** `92d7a1bef111d0dde70cd809456cd1b7`  
  **Description:** Species/group configuration for Phylo Permul Habitat Forest Savanna; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Habitat_Forest_Shrubland_config.txt`  
  **MD5:** `b1334ad2db694322015639386a94d9a9`  
  **Description:** Species/group configuration for Phylo Permul Habitat Forest Shrubland; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Habitat_Forest_Wetlands_config.txt`  
  **MD5:** `890b88721d1ffebf4a0bf373159bf10c`  
  **Description:** Species/group configuration for Phylo Permul Habitat Forest Wetlands; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Habitat_Forest_config.txt`  
  **MD5:** `6b58e8056bcc92267c386f4b66f5681d`  
  **Description:** Species/group configuration for Phylo Permul Habitat Forest; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_HeadFlesh_config.txt`  
  **MD5:** `f084c03b5e509a3626045d8114b3f37b`  
  **Description:** Species/group configuration for Phylo Permul HeadFlesh; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_HeadHair_config.txt`  
  **MD5:** `be300d642bbd120b05eda2ea8d4a53bc`  
  **Description:** Species/group configuration for Phylo Permul HeadHair; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Heart_Mass_g_config.txt`  
  **MD5:** `4d6cb8b8560e0a509b099b01f1e608ec`  
  **Description:** Species/group configuration for Phylo Permul Heart Mass g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_HomeRange_km2_config.txt`  
  **MD5:** `682dd4b4ee64eacbcc9c00b2189254aa`  
  **Description:** Species/group configuration for Phylo Permul HomeRange km2; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_IUCN_ER_config.txt`  
  **MD5:** `9fa0b87dccb03883766f2e5d79baf92d`  
  **Description:** Species/group configuration for Phylo Permul IUCN ER; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_IUCN_Pop_Trend_config.txt`  
  **MD5:** `c5812b93ef35ac34c83cc27c329b3d6a`  
  **Description:** Species/group configuration for Phylo Permul IUCN Pop Trend; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Infanticide_config.txt`  
  **MD5:** `f7097a65cb2a4763afab31e37569e38b`  
  **Description:** Species/group configuration for Phylo Permul Infanticide; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Interbirth_interval__days__config.txt`  
  **MD5:** `640347a349ad678310015cb2eabb93f4`  
  **Description:** Species/group configuration for Phylo Permul Interbirth interval days; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Lactation__days__config.txt`  
  **MD5:** `919b0eb0914efc4b26cc47455e6c431d`  
  **Description:** Species/group configuration for Phylo Permul Lactation days; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Lactic_Dehydrogenase_config.txt`  
  **MD5:** `d098a72034033c9dc601ecab938e7ca1`  
  **Description:** Species/group configuration for Phylo Permul Lactic Dehydrogenase; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Leaves_config.txt`  
  **MD5:** `70a380dc8c4c98d056a3ae8363d27e61`  
  **Description:** Species/group configuration for Phylo Permul Leaves; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_LitterSz_config.txt`  
  **MD5:** `09e43687de2a92209c716f40f1fec4a3`  
  **Description:** Species/group configuration for Phylo Permul LitterSz; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_LiverMass_g_config.txt`  
  **MD5:** `3517cc60b6bc9f152e95c2c6adb96a79`  
  **Description:** Species/group configuration for Phylo Permul LiverMass g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Locomotion_AR_config.txt`  
  **MD5:** `4fbc0b05414c71c5199a3189e4b22dcd`  
  **Description:** Species/group configuration for Phylo Permul Locomotion AR; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Locomotion_BOTH_config.txt`  
  **MD5:** `1d27bdc876083a2468fcf942a8b4340d`  
  **Description:** Species/group configuration for Phylo Permul Locomotion BOTH; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Locomotion_T_config.txt`  
  **MD5:** `8a060fe31f831f5328dbddc54bc9791c`  
  **Description:** Species/group configuration for Phylo Permul Locomotion T; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Lymphocytes_Abs_config.txt`  
  **MD5:** `1aafd62a1705873de4d1949ff5a6b787`  
  **Description:** Species/group configuration for Phylo Permul Lymphocytes Abs; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Max_All_NeonateWeight_g_config.txt`  
  **MD5:** `d50bbf3f8fb582b105e2f06802f14c59`  
  **Description:** Species/group configuration for Phylo Permul M Max All NeonateWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Max_All_YngAdultWeight_g_config.txt`  
  **MD5:** `fa3335e3f318ca5146f34d9ee20197f0`  
  **Description:** Species/group configuration for Phylo Permul M Max All YngAdultWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Max_F_YngAdultWeight_g_config.txt`  
  **MD5:** `cd4bf8001f6389e7276bfa6118238e8e`  
  **Description:** Species/group configuration for Phylo Permul M Max F YngAdultWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Max_M_YngAdultWeight_g_config.txt`  
  **MD5:** `72cf8b92ff234a189ceea11b29406191`  
  **Description:** Species/group configuration for Phylo Permul M Max M YngAdultWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Mean_M_NeonateWeight_g_config.txt`  
  **MD5:** `d4b902b89b4aa4107a52f5905458f8e8`  
  **Description:** Species/group configuration for Phylo Permul M Mean M NeonateWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_M_Min_M_NeonateWeight_g_config.txt`  
  **MD5:** `13d46eb0c3f0b5291c49552926160076`  
  **Description:** Species/group configuration for Phylo Permul M Min M NeonateWeight g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Monogamous_config.txt`  
  **MD5:** `4c7bdac0f3612684723cd185c8c96edc`  
  **Description:** Species/group configuration for Phylo Permul Mating system Monogamous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Polygynandrous_Monogamous_config.txt`  
  **MD5:** `7ab3ca8215e293bde46da0e24afbf309`  
  **Description:** Species/group configuration for Phylo Permul Mating system Polygynandrous Monogamous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Polygynandrous_config.txt`  
  **MD5:** `c346d69610c87d0cfd848aa2735f95d5`  
  **Description:** Species/group configuration for Phylo Permul Mating system Polygynandrous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Polygynous_Monogamous_config.txt`  
  **MD5:** `4ad5975ea060c325176d9c4421559408`  
  **Description:** Species/group configuration for Phylo Permul Mating system Polygynous Monogamous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Polygynous_Polygynandrous_config.txt`  
  **MD5:** `b19bd4f45e4dabea6049aa63e9e8ccdc`  
  **Description:** Species/group configuration for Phylo Permul Mating system Polygynous Polygynandrous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mating_system_Polygynous_config.txt`  
  **MD5:** `f7d79d9cd261c7d82b863cb9ef93f955`  
  **Description:** Species/group configuration for Phylo Permul Mating system Polygynous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Maximum_possible_reproductive_success_under_optimal_conditions__rmax__config.txt`  
  **MD5:** `cacf6afe1997c8da89a953efc9993732`  
  **Description:** Species/group configuration for Phylo Permul Maximum possible reproductive success under optimal conditions rmax; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_MeanGroupSize_config.txt`  
  **MD5:** `c2bfb11117bee9b244a049723f240e35`  
  **Description:** Species/group configuration for Phylo Permul MeanGroupSize; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Mesencephalon_config.txt`  
  **MD5:** `0dc9e2d347ad8fa0097a6d20f81aee4d`  
  **Description:** Species/group configuration for Phylo Permul Mesencephalon; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Monocytes_Abs_config.txt`  
  **MD5:** `caf9ab845f9e31ca986e78f38f92aed5`  
  **Description:** Species/group configuration for Phylo Permul Monocytes Abs; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Monocytes___config.txt`  
  **MD5:** `14e3ac1ad7d1d13ad5921ff2d4a777ba`  
  **Description:** Species/group configuration for Phylo Permul Monocytes; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Monopolization_Yes_config.txt`  
  **MD5:** `b8421f5042e853e6a0a33b20403b560c`  
  **Description:** Species/group configuration for Phylo Permul Monopolization Yes; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Move_config.txt`  
  **MD5:** `8d37482dc3d1094d0bccb15d82e9b66a`  
  **Description:** Species/group configuration for Phylo Permul Move; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_NatalDispersal_Sexbiased_config.txt`  
  **MD5:** `8b077a30473c410c6be2a3acd7c098d2`  
  **Description:** Species/group configuration for Phylo Permul NatalDispersal Sexbiased; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Neutrophil_Band_Abs_config.txt`  
  **MD5:** `75de35b16235f7807a33d639c854f7a7`  
  **Description:** Species/group configuration for Phylo Permul Neutrophil Band Abs; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Neutrophil_Band___config.txt`  
  **MD5:** `75de35b16235f7807a33d639c854f7a7`  
  **Description:** Species/group configuration for Phylo Permul Neutrophil Band; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_PET_Mean_mm_config.txt`  
  **MD5:** `b1d7b7897f376f085317b7a0d999c6ee`  
  **Description:** Species/group configuration for Phylo Permul PET Mean mm; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Paternal_care_config.txt`  
  **MD5:** `aee839af1a3c2a5848b03f7d96801603`  
  **Description:** Species/group configuration for Phylo Permul Paternal care; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Positional_behaviour_config.txt`  
  **MD5:** `d18381da6d735a595391e30e99d6f827`  
  **Description:** Species/group configuration for Phylo Permul Positional behaviour; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Precip_Mean_mm_config.txt`  
  **MD5:** `b57b1b876b80eab2d2de71fe8521dd83`  
  **Description:** Species/group configuration for Phylo Permul Precip Mean mm; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_R_Max_LitterSize_config.txt`  
  **MD5:** `292d48ca1609853e77d8ecea44155a71`  
  **Description:** Species/group configuration for Phylo Permul R Max LitterSize; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_R_Min_Dam_AgeAtConcep_y_config.txt`  
  **MD5:** `cf4d583dcb1239e7ff8265663f6b7f97`  
  **Description:** Species/group configuration for Phylo Permul R Min Dam AgeAtConcep y; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Seasonal_breeders_config.txt`  
  **MD5:** `ea2507c33c51b388b2e8715d2fe10aba`  
  **Description:** Species/group configuration for Phylo Permul Seasonal breeders; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Septum_config.txt`  
  **MD5:** `ffe9c6343f4ca573f58a41bdf9fd2fac`  
  **Description:** Species/group configuration for Phylo Permul Septum; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Ser_Glut_Pyruvic_Transaminase_config.txt`  
  **MD5:** `3c2c976a28b14ba2cbfc7f822f4692aa`  
  **Description:** Species/group configuration for Phylo Permul Ser Glut Pyruvic Transaminase; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SexualOrnamentDimorphism__ExcludingDichromatism__config.txt`  
  **MD5:** `dc083b766e6b741882d8d645d14f93d9`  
  **Description:** Species/group configuration for Phylo Permul SexualOrnamentDimorphism ExcludingDichromatism; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SexualOrnamentDimorphism_config.txt`  
  **MD5:** `b67d40f1067abbba3df810f71461dde5`  
  **Description:** Species/group configuration for Phylo Permul SexualOrnamentDimorphism; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SocialMatingSystem_Monogamous_config.txt`  
  **MD5:** `91c9b295d1a9bb4069ed2c9a6db80a2b`  
  **Description:** Species/group configuration for Phylo Permul SocialMatingSystem Monogamous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SocialMatingSystem_Multimale_config.txt`  
  **MD5:** `ecb1ffd78160a0f8b05f5db5a14369ab`  
  **Description:** Species/group configuration for Phylo Permul SocialMatingSystem Multimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SocialMatingSystem_Polygynous_config.txt`  
  **MD5:** `b8421f5042e853e6a0a33b20403b560c`  
  **Description:** Species/group configuration for Phylo Permul SocialMatingSystem Polygynous; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_config.txt`  
  **MD5:** `8a2d38f066fa57123b9708a29119bed8`  
  **Description:** Species/group configuration for Phylo Permul Social; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Multimale_config.txt`  
  **MD5:** `dd94c1216c4e2d06861c22870a427612`  
  **Description:** Species/group configuration for Phylo Permul Social organization Multimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Pair_Multimale_config.txt`  
  **MD5:** `14ee10791258d959f5bbd32328205fd3`  
  **Description:** Species/group configuration for Phylo Permul Social organization Pair Multimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Pair_Unimale_Multimale_config.txt`  
  **MD5:** `04a63e26a06653fb75f7c8047a9d9113`  
  **Description:** Species/group configuration for Phylo Permul Social organization Pair Unimale Multimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Pair_Unimale_config.txt`  
  **MD5:** `1835b83974c6a9b3129460537e15a073`  
  **Description:** Species/group configuration for Phylo Permul Social organization Pair Unimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Pair_config.txt`  
  **MD5:** `e60dd4f5b45dba14874f13391e14e8b2`  
  **Description:** Species/group configuration for Phylo Permul Social organization Pair; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Solitary_config.txt`  
  **MD5:** `e161c4da7742a9f6cc491669124ab590`  
  **Description:** Species/group configuration for Phylo Permul Social organization Solitary; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Unimale_Multimale_config.txt`  
  **MD5:** `9c25c76177fba5b42de8b8c3df79c7a1`  
  **Description:** Species/group configuration for Phylo Permul Social organization Unimale Multimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Social_organization_Unimale_config.txt`  
  **MD5:** `a3780976315f476e289c509d062c02cc`  
  **Description:** Species/group configuration for Phylo Permul Social organization Unimale; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Sodium_config.txt`  
  **MD5:** `8f2be7f72fe29d8b00c35819e06f9b3c`  
  **Description:** Species/group configuration for Phylo Permul Sodium; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_SpleenMass_g_config.txt`  
  **MD5:** `73b8f967b222cd00d1664f240dd3e91d`  
  **Description:** Species/group configuration for Phylo Permul SpleenMass g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_StableGroups_config.txt`  
  **MD5:** `3dda6fbc8a589b6e7e653c6f58070cb9`  
  **Description:** Species/group configuration for Phylo Permul StableGroups; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_StomachMass_g_config.txt`  
  **MD5:** `743c5d0920aeb8e6feceabc551a48183`  
  **Description:** Species/group configuration for Phylo Permul StomachMass g; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Territoriality_config.txt`  
  **MD5:** `e92b9089d772a0c6684b909e18a7aea2`  
  **Description:** Species/group configuration for Phylo Permul Territoriality; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TotalOrnamentDimorphism_config.txt`  
  **MD5:** `0b1af4b8eaae1b29841389b195c2b383`  
  **Description:** Species/group configuration for Phylo Permul TotalOrnamentDimorphism; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_Triglycerides_config.txt`  
  **MD5:** `77c46a09d3b5907b06906a25e8742cd8`  
  **Description:** Species/group configuration for Phylo Permul Triglycerides; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Folivore_config.txt`  
  **MD5:** `34bbfd1211a7cb6dd887e8af70c873b6`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Folivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Folivore_frugivore_config.txt`  
  **MD5:** `5f0fa62b9da1bd79f04af04172e469b5`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Folivore frugivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Frugivore_config.txt`  
  **MD5:** `b433687f86b74388de71bb5f53a8a02f`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Frugivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Gummivore_config.txt`  
  **MD5:** `530ca6b79c7dadd852460a73c839c0da`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Gummivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Insectivore_config.txt`  
  **MD5:** `cd188f92915a9d680559a3ee7c84285a`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Insectivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_TrophicGuild_Omnivore_config.txt`  
  **MD5:** `0be0f01eea4b9e5882fcf7be9665a060`  
  **Description:** Species/group configuration for Phylo Permul TrophicGuild Omnivore; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Phylo_Permul_White_Blood_Cells_config.txt`  
  **MD5:** `2aec7f968dd5476620ea72d3ff714bf1`  
  **Description:** Species/group configuration for Phylo Permul White Blood Cells; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Piriform_Lobe_config.txt`  
  **MD5:** `37ed195894569db8258861afe4ac7762`  
  **Description:** Species/group configuration for Piriform Lobe; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Schizocortex_config.txt`  
  **MD5:** `70334b42efa0a777331ca1c43361f4b8`  
  **Description:** Species/group configuration for Schizocortex; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Ser_Glut_Oxaloacetic_Transamin__config.txt`  
  **MD5:** `f3353d5d133761c2ee34727e06ea8faa`  
  **Description:** Species/group configuration for Ser Glut Oxaloacetic Transamin; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Striatum__incl__NAcc__config.txt`  
  **MD5:** `08ce225984fef58e30190f67c7e4aefd`  
  **Description:** Species/group configuration for Striatum incl NAcc; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Striatum_config.txt`  
  **MD5:** `ab8055bbe11902444e9d44b99ce76fbb`  
  **Description:** Species/group configuration for Striatum; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Thalamus_config.txt`  
  **MD5:** `b1a57987141fcce4a8b118cc07a6ea96`  
  **Description:** Species/group configuration for Thalamus; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/Uric_Acid_config.txt`  
  **MD5:** `dfd75f9c0317862387b80eb9ad114b04`  
  **Description:** Species/group configuration for Uric Acid; used to define or compare trait groups.

- **File:** `workflow/SUPP/S27/Configs_raw/XII_config.txt`  
  **MD5:** `a6a330648ba1f4b06bfe99ca42dfa95e`  
  **Description:** Species/group configuration for XII; used to define or compare trait groups.

### workflow/SUPP/S27

- **File:** `workflow/SUPP/S27/FigS27.R`  
  **MD5:** `bec4079cb44aa01316d38f98b29c3f10`  
  **Description:** R analysis and plotting script for Supplementary Figure 27.

### workflow/SUPP/S28

- **File:** `workflow/SUPP/S28/ALL_PHYLOGISTIC_OCT23_BETAS.txt`  
  **MD5:** `757fda0c13106901d1d34b5ecb85a329`  
  **Description:** Phylogenetic logistic model association results.

- **File:** `workflow/SUPP/S28/FigS28.R`  
  **MD5:** `81dae8179b1c5674562368e28809ca81`  
  **Description:** R analysis and plotting script for Supplementary Figure 28.

### workflow/SUPP/S29

- **File:** `workflow/SUPP/S29/ALL_POWER_TRAITS_MAR23.txt`  
  **MD5:** `bc6095c99a32da1871d01aa122be6ea0`  
  **Description:** Trait-level association/power-analysis input used by supplementary scripts.

- **File:** `workflow/SUPP/S29/FigS29.R`  
  **MD5:** `942bedb06b8ca1cf34941787c237bcff`  
  **Description:** R analysis and plotting script for Supplementary Figure 29.

- **File:** `workflow/SUPP/S29/Supp2Table3.tsv`  
  **MD5:** `1be14af08da78c5eb09f15a3d2d5d8c0`  
  **Description:** Trait contrast, family, outlier, species and permutation-count metadata.

### workflow/SUPP/S3

- **File:** `workflow/SUPP/S3/FigS3.R`  
  **MD5:** `ff5c4cc75165eac030a86bb6222fea25`  
  **Description:** R analysis and plotting script for Supplementary Figure 3.

### workflow/SUPP/S30

- **File:** `workflow/SUPP/S30/FigS30.R`  
  **MD5:** `f72ef4bbe8ce0d6edc9cceda00e2bf89`  
  **Description:** R analysis and plotting script for Supplementary Figure 30.

### workflow/SUPP/S31

- **File:** `workflow/SUPP/S31/FigS31.R`  
  **MD5:** `ab660e1d58cc64cbc3e3a80875283783`  
  **Description:** R analysis and plotting script for Supplementary Figure 31.

### workflow/SUPP/S32

- **File:** `workflow/SUPP/S32/AAs_colors.gz`  
  **MD5:** `141f8c5da8d5c18d697131abda5d722f`  
  **Description:** Compressed serialized R amino-acid color palette.

- **File:** `workflow/SUPP/S32/BLOSUM62.txt`  
  **MD5:** `7e0c9e55df34f8002e9a1ae6cfaf087f`  
  **Description:** BLOSUM62 amino-acid substitution matrix.

- **File:** `workflow/SUPP/S32/FigS32.R`  
  **MD5:** `95725bfe1c147fd165475546c9aaac95`  
  **Description:** R analysis and plotting script for Supplementary Figure 32.

- **File:** `workflow/SUPP/S32/Table_AAs_TOP.txt`  
  **MD5:** `240383b92c2481f85ac8df3e4c6c3879`  
  **Description:** Amino-acid table used by Supplementary Figure 32.

### workflow/SUPP/S33

- **File:** `workflow/SUPP/S33/FigS33.R`  
  **MD5:** `d0c73fd91b609490b6bd0491e7e9067b`  
  **Description:** R analysis and plotting script for Supplementary Figure 33.

- **File:** `workflow/SUPP/S33/Table_CAAs_SCENARIO1_2.txt`  
  **MD5:** `2a1df0083d5ab270108359b7d208b9c7`  
  **Description:** CAAS entries labeled by scenario, including gene, trait and residue information.

### workflow/SUPP/S34

- **File:** `workflow/SUPP/S34/ALL_MAMMALIAN_INFO_MAY23_FILTTOT.txt`  
  **MD5:** `b1ca56e9dba8847722ea0d433613e102`  
  **Description:** Filtered mammalian positional and annotation information.

- **File:** `workflow/SUPP/S34/ALL_TRANSVAR_INFO_APR23.txt`  
  **MD5:** `6a296ddbb6ba37c9f1a929b54df227ab`  
  **Description:** Gene, trait, position and TransVar mutation annotations.

- **File:** `workflow/SUPP/S34/FigS34.R`  
  **MD5:** `879021609651253a2cd2bc9f4e53080b`  
  **Description:** R analysis and plotting script for Supplementary Figure 34.

### workflow/SUPP/S35

- **File:** `workflow/SUPP/S35/FigS35.R`  
  **MD5:** `04cba3e137a253b474082541d0fec0b6`  
  **Description:** R analysis and plotting script for Supplementary Figure 35.

### workflow/SUPP/S36

- **File:** `workflow/SUPP/S36/FigS36.R`  
  **MD5:** `b3ed3bc33d6e2c982388f8fcc7a47b01`  
  **Description:** R analysis and plotting script for Supplementary Figure 36.

### workflow/SUPP/S36/PHYLOP_files

- **File:** `workflow/SUPP/S36/PHYLOP_files/ALL_JOINED.gz`  
  **MD5:** `ed7b300460a967cdb4227e6d3faf0a32`  
  **Description:** Compressed serialized R phyloP analysis object: ALL JOINED. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/ALL_JOINED_CONTR.gz`  
  **MD5:** `4598ea7ca49a581b0a593c8539eef628`  
  **Description:** Compressed serialized R phyloP analysis object: ALL JOINED CONTR. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/SCENARIO1.gz`  
  **MD5:** `4242bb9f23f0e1678f1abca6dd7efe34`  
  **Description:** Compressed serialized R phyloP analysis object: SCENARIO1. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/SCENARIO1_CONTR.gz`  
  **MD5:** `39048bc3db04b6d30e3fdf3373d66639`  
  **Description:** Compressed serialized R phyloP analysis object: SCENARIO1 CONTR. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/SCENARIO2.gz`  
  **MD5:** `3c29311843892fab2e49870b7ddcc591`  
  **Description:** Compressed serialized R phyloP analysis object: SCENARIO2. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/SCENARIO2_CONTR.gz`  
  **MD5:** `3807b62bcac74eecc48942ebefa65f69`  
  **Description:** Compressed serialized R phyloP analysis object: SCENARIO2 CONTR. Scenario labels are retained from the source.

- **File:** `workflow/SUPP/S36/PHYLOP_files/Sim.gz`  
  **MD5:** `9920b99ff6d54d19ea5c291804f3214e`  
  **Description:** Compressed serialized R phyloP analysis object: Sim. Scenario labels are retained from the source.

### workflow/SUPP/S36/PHYLOP_files/all_phyloP

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/ALL_CDS_PHYLOP_shuff1.out`  
  **MD5:** `02095138aed7f4927abb3e56c6eac95c`  
  **Description:** PhyloP output for shuffled coding-sequence positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr1.phyloPs.out`  
  **MD5:** `bfd13da92d56b3f6777a8866ee649b5a`  
  **Description:** Chromosome 1 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr10.phyloPs.out`  
  **MD5:** `a76a18d2fb763108916cb0b2b579cd93`  
  **Description:** Chromosome 10 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr13.phyloPs.out`  
  **MD5:** `c8b524954355d0d8f2bf5cff65ace689`  
  **Description:** Chromosome 13 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr14.phyloPs.out`  
  **MD5:** `c8a69e13d7b9502f60cd08751841d12a`  
  **Description:** Chromosome 14 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr15.phyloPs.out`  
  **MD5:** `000a7c7048efb974c8fee1facbc97c9b`  
  **Description:** Chromosome 15 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr17.phyloPs.out`  
  **MD5:** `489bd7b56a770481f8f4aaba3bf3ee68`  
  **Description:** Chromosome 17 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr18.phyloPs.out`  
  **MD5:** `db9365ddf12bd3ac59e5538581d88c5d`  
  **Description:** Chromosome 18 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr19.phyloPs.out`  
  **MD5:** `80a2ed20bec2b5d0dbd52a1b44eedc62`  
  **Description:** Chromosome 19 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr2.phyloPs.out`  
  **MD5:** `5e3dae524237cb98833498b5d8fcd391`  
  **Description:** Chromosome 2 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr21.phyloPs.out`  
  **MD5:** `5af504214fbcd7c854654c3dba67360f`  
  **Description:** Chromosome 21 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr22.phyloPs.out`  
  **MD5:** `ecba802e6a81b4efec28f6d8080a5975`  
  **Description:** Chromosome 22 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr3.phyloPs.out`  
  **MD5:** `cbd2ed03c0bb8041d7fa07290c0e506e`  
  **Description:** Chromosome 3 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr4.phyloPs.out`  
  **MD5:** `d9585a0eb2c44b09475fa86f7fe06d0b`  
  **Description:** Chromosome 4 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr5.phyloPs.out`  
  **MD5:** `db69d0383190b4ec23974d7fb64a3be2`  
  **Description:** Chromosome 5 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr6.phyloPs.out`  
  **MD5:** `9a8780bfccf6242e8ae8b5305339f2b2`  
  **Description:** Chromosome 6 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr8.phyloPs.out`  
  **MD5:** `e3e2c7282939806132a220aa010ea5c7`  
  **Description:** Chromosome 8 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chr9.phyloPs.out`  
  **MD5:** `448bc49e692b33cc97c6d19ba829eaf2`  
  **Description:** Chromosome 9 phyloP output for comparison positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/all_phyloP/chrX.phyloPs.out`  
  **MD5:** `16f1a7bcd41f217f73240621c5a06d0b`  
  **Description:** Chromosome X phyloP output for comparison positions.

### workflow/SUPP/S36/PHYLOP_files/caas_phyloP

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr10.phyloPs.out`  
  **MD5:** `288de7f6a5b743d2a3b49f15f2d9e478`  
  **Description:** Chromosome 10 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr13.phyloPs.out`  
  **MD5:** `bd0111531d14a479a21bb33d3a0a21b8`  
  **Description:** Chromosome 13 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr14.phyloPs.out`  
  **MD5:** `11f70d09249f1fb46f2e5e4684760fc4`  
  **Description:** Chromosome 14 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr15.phyloPs.out`  
  **MD5:** `5f2d723be6042a221c98b9b4efaf4c45`  
  **Description:** Chromosome 15 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr18.phyloPs.out`  
  **MD5:** `64386fe516ecab3916b7bd721ad84e3a`  
  **Description:** Chromosome 18 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr20.phyloPs.out`  
  **MD5:** `1104df2ca71d6dd307e1989e9868afb9`  
  **Description:** Chromosome 20 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr21.phyloPs.out`  
  **MD5:** `8abdf0d355be36e25bfaf78f0acaffc0`  
  **Description:** Chromosome 21 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr22.phyloPs.out`  
  **MD5:** `a0d566c38898226e33dd5ba5f6fc9df0`  
  **Description:** Chromosome 22 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr4.phyloPs.out`  
  **MD5:** `a830048eb6347f2f6fc668b7f044f84c`  
  **Description:** Chromosome 4 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr7.phyloPs.out`  
  **MD5:** `10500f40c4d5ec4c0de3936a3aed0427`  
  **Description:** Chromosome 7 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr8.phyloPs.out`  
  **MD5:** `c8a4b40c9e69063b8c63ca1fd2a91d4b`  
  **Description:** Chromosome 8 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chr9.phyloPs.out`  
  **MD5:** `b06b9895e5b514bd7c054aed86a1098a`  
  **Description:** Chromosome 9 phyloP output for CAAS positions.

- **File:** `workflow/SUPP/S36/PHYLOP_files/caas_phyloP/chrX.phyloPs.out`  
  **MD5:** `15e14afc36c30113365289c8d9cc8105`  
  **Description:** Chromosome X phyloP output for CAAS positions.

### workflow/SUPP/S37

- **File:** `workflow/SUPP/S37/ALL_TRANSVAR_RESULTS_APR23_FORMATTED_FILT.txt`  
  **MD5:** `8db4f242d9733049fde9cd0fe6e642a5`  
  **Description:** Formatted and filtered genomic TransVar annotations.

- **File:** `workflow/SUPP/S37/FigS37.R`  
  **MD5:** `c5976952a6f3de95fa8c7d34210b6a16`  
  **Description:** R analysis and plotting script for Supplementary Figure 37.

### workflow/SUPP/S37/VARIATION_files

- **File:** `workflow/SUPP/S37/VARIATION_files/ALL_VARIATION_NOV22.txt`  
  **MD5:** `8bf60cac42a92a5c0ba03bb905394734`  
  **Description:** Combined variation records labeled NOV22 in the source filename.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr1.variation_APR23.txt`  
  **MD5:** `9fe439d9ea6ea669cce243a953c9a308`  
  **Description:** Chromosome 1 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr1.variation_NOV22.txt`  
  **MD5:** `aefe1f55b620e06a657d7559d9af2803`  
  **Description:** Chromosome 1 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr10.variation_APR23.txt`  
  **MD5:** `00dea4f47c1680ba0e7298330b52163b`  
  **Description:** Chromosome 10 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr10.variation_NOV22.txt`  
  **MD5:** `a6f328246aa86c74023f5dd3e37de75c`  
  **Description:** Chromosome 10 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr11.variation_APR23.txt`  
  **MD5:** `bb7a93f3971fccb541caf4b483fc352d`  
  **Description:** Chromosome 11 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr11.variation_NOV22.txt`  
  **MD5:** `b5c13abe77daa70f6336fba4b204acd3`  
  **Description:** Chromosome 11 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr12.variation_APR23.txt`  
  **MD5:** `db370e399262da149b4bd6ab944121aa`  
  **Description:** Chromosome 12 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr12.variation_NOV22.txt`  
  **MD5:** `02371d429ae720219ccb57473daca455`  
  **Description:** Chromosome 12 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr13.variation_APR23.txt`  
  **MD5:** `c9729a1d4b04641fa377daa0e7e2d89e`  
  **Description:** Chromosome 13 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr13.variation_NOV22.txt`  
  **MD5:** `a2b2894f650c16b3d8c812358448a91b`  
  **Description:** Chromosome 13 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr14.variation_APR23.txt`  
  **MD5:** `1e653ddbcc5c35cacf0a83e205d9d12c`  
  **Description:** Chromosome 14 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr14.variation_NOV22.txt`  
  **MD5:** `2d22b66754dbdbf6041e43dde299fac5`  
  **Description:** Chromosome 14 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr15.variation_APR23.txt`  
  **MD5:** `dd250e1a1cfee350053ffc436a08a029`  
  **Description:** Chromosome 15 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr15.variation_NOV22.txt`  
  **MD5:** `b40f881ccc8af0c131a64ed76bb97f27`  
  **Description:** Chromosome 15 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr16.variation_APR23.txt`  
  **MD5:** `100625045a03e53440cfba598f7dc103`  
  **Description:** Chromosome 16 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr16.variation_NOV22.txt`  
  **MD5:** `298f298900c6f12f4521352a283cb53e`  
  **Description:** Chromosome 16 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr17.variation_APR23.txt`  
  **MD5:** `c13e23cbbf5d97492f1d2cbaa36b233c`  
  **Description:** Chromosome 17 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr17.variation_NOV22.txt`  
  **MD5:** `10d7953beda7a8a31c9fadcfd5717e33`  
  **Description:** Chromosome 17 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr18.variation_APR23.txt`  
  **MD5:** `5b1c7d0695d40585996ff6dce26b3749`  
  **Description:** Chromosome 18 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr18.variation_NOV22.txt`  
  **MD5:** `37c2978614a2039005195e3a1ee440c6`  
  **Description:** Chromosome 18 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr19.variation_APR23.txt`  
  **MD5:** `f9de54ff61eceb04d20a2be6a770a70b`  
  **Description:** Chromosome 19 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr19.variation_NOV22.txt`  
  **MD5:** `91659733ac145e7231212317d1f1b486`  
  **Description:** Chromosome 19 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr2.variation_APR23.txt`  
  **MD5:** `5c73cb489b906f4ff937fa2b99901282`  
  **Description:** Chromosome 2 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr2.variation_NOV22.txt`  
  **MD5:** `4b80f71a41363d8973e66fde739e61ed`  
  **Description:** Chromosome 2 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr20.variation_APR23.txt`  
  **MD5:** `89443ef43c8afa2c9fdf73f9c659a6b5`  
  **Description:** Chromosome 20 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr20.variation_NOV22.txt`  
  **MD5:** `94e4e0a7a40c2d961ff5d8862fdb07b5`  
  **Description:** Chromosome 20 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr21.variation_APR23.txt`  
  **MD5:** `7e7dad686f700c5e9eb785c4b2cbd430`  
  **Description:** Chromosome 21 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr21.variation_NOV22.txt`  
  **MD5:** `9bf92865f6465aef13b2c185e50a2e28`  
  **Description:** Chromosome 21 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr22.variation_APR23.txt`  
  **MD5:** `6f67f42ccb55fa283930e067d0436d6a`  
  **Description:** Chromosome 22 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr22.variation_NOV22.txt`  
  **MD5:** `01cc83a42dd874347d09517079f2c2d6`  
  **Description:** Chromosome 22 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr3.variation_APR23.txt`  
  **MD5:** `1c920afcd13e033bffb0e9e483b8532d`  
  **Description:** Chromosome 3 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr3.variation_NOV22.txt`  
  **MD5:** `a6af736479433c77684f38a9136914dc`  
  **Description:** Chromosome 3 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr4.variation_APR23.txt`  
  **MD5:** `dd9aa34932562f1e592d31b744e32584`  
  **Description:** Chromosome 4 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr4.variation_NOV22.txt`  
  **MD5:** `2de011a3e43049c9485755dda69fff12`  
  **Description:** Chromosome 4 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr5.variation_APR23.txt`  
  **MD5:** `e9a93369c6e23c6e5bd5562935b04bda`  
  **Description:** Chromosome 5 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr5.variation_NOV22.txt`  
  **MD5:** `822c9079800f1c0e8221119f145db645`  
  **Description:** Chromosome 5 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr6.variation_APR23.txt`  
  **MD5:** `3bae7551079af32ae4dcb949dce5f172`  
  **Description:** Chromosome 6 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr6.variation_NOV22.txt`  
  **MD5:** `3d30a07d5da61078b0e851320352f3ab`  
  **Description:** Chromosome 6 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr7.variation_APR23.txt`  
  **MD5:** `ebc57eb94024a1c9557bff3ab0cbdceb`  
  **Description:** Chromosome 7 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr7.variation_NOV22.txt`  
  **MD5:** `0da76f089d0ac7a21749b61a6d624a50`  
  **Description:** Chromosome 7 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr8.variation_APR23.txt`  
  **MD5:** `4338cc0d8717e2e3d1417f322c0706c2`  
  **Description:** Chromosome 8 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr8.variation_NOV22.txt`  
  **MD5:** `cf346609634896a58930b19ff2baec47`  
  **Description:** Chromosome 8 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr9.variation_APR23.txt`  
  **MD5:** `cf4b96fff22f3e5992e049f977e154e8`  
  **Description:** Chromosome 9 variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chr9.variation_NOV22.txt`  
  **MD5:** `6ca0c2e98c0ba53a40ad561130434fac`  
  **Description:** Chromosome 9 variation records; source version NOV22.

- **File:** `workflow/SUPP/S37/VARIATION_files/chrX.variation_APR23.txt`  
  **MD5:** `623a631b69093c5e53cba18ede5af50d`  
  **Description:** Chromosome X variation records; source version APR23.

- **File:** `workflow/SUPP/S37/VARIATION_files/chrX.variation_NOV22.txt`  
  **MD5:** `c3ec5a3324c99a2ce33d9aba7e9a2857`  
  **Description:** Chromosome X variation records; source version NOV22.

### workflow/SUPP/S38

- **File:** `workflow/SUPP/S38/FigS38.R`  
  **MD5:** `a7e3e7d93df79a3df6dd929aa5eea593`  
  **Description:** R analysis and plotting script for Supplementary Figure 38.

- **File:** `workflow/SUPP/S38/joined_FEB24_QUAL.gz`  
  **MD5:** `7223cc4b27cc4eadcc7c0c0b70d182af`  
  **Description:** Compressed serialized R joined qualitative-trait results for Supplementary Figure 38.

### workflow/SUPP/S39

- **File:** `workflow/SUPP/S39/FigS39.R`  
  **MD5:** `db8d7820fa3f1fbc2c2ba113cf064826`  
  **Description:** R analysis and plotting script for Supplementary Figure 39.

### workflow/SUPP/S4

- **File:** `workflow/SUPP/S4/FigS4.R`  
  **MD5:** `978215a88e731ed377db7a23cff2fcfe`  
  **Description:** R analysis and plotting script for Supplementary Figure 4.

- **File:** `workflow/SUPP/S4/Redundant_measurements.tsv`  
  **MD5:** `32897da2b11f9d59431be7a6564b9779`  
  **Description:** Trait measurement redundancy list.

### workflow/SUPP/S40

- **File:** `workflow/SUPP/S40/ALL_POWER_TRAITS_MAR23.txt`  
  **MD5:** `bc6095c99a32da1871d01aa122be6ea0`  
  **Description:** Trait-level association/power-analysis input used by supplementary scripts.

- **File:** `workflow/SUPP/S40/FigS40.R`  
  **MD5:** `f1610f9207aff3ace2fc9067384cff93`  
  **Description:** R analysis and plotting script for Supplementary Figure 40.

- **File:** `workflow/SUPP/S40/SPECIES_COUNTS_TRAITS_PASTED.txt`  
  **MD5:** `c3ec54127d91c729f1bb002b3f4478e2`  
  **Description:** Species-count metadata by trait used in power-analysis figures.

### workflow/SUPP/S41

- **File:** `workflow/SUPP/S41/FigS41.R`  
  **MD5:** `ade80c13737aa230a890a08d58cc4671`  
  **Description:** R analysis and plotting script for Supplementary Figure 41.

### workflow/SUPP/S42

- **File:** `workflow/SUPP/S42/FigS42.R`  
  **MD5:** `4c58ff23fcc954e76d11504130fb2dea`  
  **Description:** R analysis and plotting script for Supplementary Figure 42.

### workflow/SUPP/S43

- **File:** `workflow/SUPP/S43/FigS43.R`  
  **MD5:** `e1991093f0969a0af113d918f34eaaf2`  
  **Description:** R analysis and plotting script for Supplementary Figure 43.

### workflow/SUPP/S44

- **File:** `workflow/SUPP/S44/FigS44.R`  
  **MD5:** `541b0dbe7b1b09ae491bd4b96bf980bc`  
  **Description:** R analysis and plotting script for Supplementary Figure 44.

### workflow/SUPP/S45

- **File:** `workflow/SUPP/S45/ALL_BM_APR23.txt`  
  **MD5:** `d013f4908bd8d26f9ccc2f6680f300a1`  
  **Description:** Brownian-motion association results used by Supplementary Figure 45.

- **File:** `workflow/SUPP/S45/ALL_LAMBDA_APR23.txt`  
  **MD5:** `a40a005b99f03140be0142065d0e71f5`  
  **Description:** Lambda model association results used by Supplementary Figure 45.

- **File:** `workflow/SUPP/S45/ALL_POWER_APR23.txt`  
  **MD5:** `a6986d061daad16c883edfcb561513de`  
  **Description:** Power-analysis results used by Supplementary Figure 45.

- **File:** `workflow/SUPP/S45/FigS45.R`  
  **MD5:** `95b2a7535e2dbd87a4bd594c62671d9a`  
  **Description:** R analysis and plotting script for Supplementary Figure 45.

### workflow/SUPP/S46

- **File:** `workflow/SUPP/S46/FigS46.R`  
  **MD5:** `8cd5f4178240220f7c5a8e7edb1d4785`  
  **Description:** R analysis and plotting script for Supplementary Figure 46.

### workflow/SUPP/S47

- **File:** `workflow/SUPP/S47/FigS47.R`  
  **MD5:** `96bd245a774f27b31ca8e7b779fee88a`  
  **Description:** R analysis and plotting script for Supplementary Figure 47.

### workflow/SUPP/S48

- **File:** `workflow/SUPP/S48/FigS48.R`  
  **MD5:** `5166a1076195d63e039960979ab8658e`  
  **Description:** R analysis and plotting script for Supplementary Figure 48.

### workflow/SUPP/S49

- **File:** `workflow/SUPP/S49/FigS49.R`  
  **MD5:** `30e9e079de40df6775cce5a36d2ad221`  
  **Description:** R analysis and plotting script for Supplementary Figure 49.

### workflow/SUPP/S5

- **File:** `workflow/SUPP/S5/FigS5.R`  
  **MD5:** `314c59f9c686569ff68be10c1ae3d6eb`  
  **Description:** R analysis and plotting script for Supplementary Figure 5.

### workflow/SUPP/S50

- **File:** `workflow/SUPP/S50/FigS50.R`  
  **MD5:** `f0beeeb0bb561e368d959f62896daaae`  
  **Description:** R analysis and plotting script for Supplementary Figure 50.

### workflow/SUPP/S6

- **File:** `workflow/SUPP/S6/FigS6.R`  
  **MD5:** `40ab5686db68d52c6d11d8e9499479ce`  
  **Description:** R analysis and plotting script for Supplementary Figure 6.

### workflow/SUPP/S6,9,12,15,18,21_MERGED

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S12_MERGED.R`  
  **MD5:** `5afd176caf1b1f789ff5244ee3c72d94`  
  **Description:** R analysis and plotting script for Supplementary Figure 12; merged analysis variant; retained syntax error.

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S15_MERGED.R`  
  **MD5:** `c3cb8e575e8f06230a3f951c8dc5e56c`  
  **Description:** R analysis and plotting script for Supplementary Figure 15; merged analysis variant; retained syntax error.

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S18_MERGED.R`  
  **MD5:** `24a516e34c6a774168ef6f315ba16798`  
  **Description:** R analysis and plotting script for Supplementary Figure 18; merged analysis variant; retained syntax error.

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S21_MERGED.R`  
  **MD5:** `f9e3be9ae1d27320c49e4207f6e1ed84`  
  **Description:** R analysis and plotting script for Supplementary Figure 21; merged analysis variant.

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S6_MERGED.R`  
  **MD5:** `707ac15bb8d1409a25628857311e7d7b`  
  **Description:** R analysis and plotting script for Supplementary Figure 6; merged analysis variant; retained syntax error.

- **File:** `workflow/SUPP/S6,9,12,15,18,21_MERGED/S9_MERGED.R`  
  **MD5:** `08a649d6c740b35c7bc9617fd99815d0`  
  **Description:** R analysis and plotting script for Supplementary Figure 9; merged analysis variant; retained syntax error.

### workflow/SUPP/S7

- **File:** `workflow/SUPP/S7/FigS7.R`  
  **MD5:** `4b715c49550f054d3c79abf00eaa901f`  
  **Description:** R analysis and plotting script for Supplementary Figure 7.

### workflow/SUPP/S8

- **File:** `workflow/SUPP/S8/FigS8.R`  
  **MD5:** `946d7d052c2a1c6c0d64637a39c270b9`  
  **Description:** R analysis and plotting script for Supplementary Figure 8.

### workflow/SUPP/S9

- **File:** `workflow/SUPP/S9/FigS9.R`  
  **MD5:** `6a712899dc01f54769f59299f56a6fc5`  
  **Description:** R analysis and plotting script for Supplementary Figure 9.

### workflow/SUPP/s11

- **File:** `workflow/SUPP/s11/FigS11.R`  
  **MD5:** `41bea09bd73838a0288e8aeabab009c7`  
  **Description:** R analysis and plotting script for Supplementary Figure 11.
