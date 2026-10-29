# The Primate Genome-Phenome Archive (PGA)

The **Primate Genome-Phenome Archive (PGA)** brings together phenotypic observations, trait metadata, protein-coding sequence alignments, genome–phenome association results, and analysis materials for comparative research in primates. The collection supports the study of how molecular evolution relates to variation in organismal traits, including analyses of convergent amino acid substitutions and relative evolutionary rates.

This repository accompanies *Valenzuela et al. (2026), A phylogenetic protein-coding genome-phenome map of complex traits across 224 primate species*. It provides access to the data used to describe phenotypic variation, define trait contrasts, investigate genetic associations, and prepare the study's figures and supplementary datasets.

The materials are organized into four main folders. `phenotype.tables/` contains the phenotypic resource and its metadata; `genome.phenome.maps/` contains the principal molecular association tables; `orthologous.alignments/` documents access to the versioned alignment collection; and `valenzuela2026.supplementary.information/` brings together supplementary workbooks, figure-analysis materials, and the archived genome–phenome workflows.

## Repository Structure

The tree below shows the main directories and selected nested components. Figure-specific inputs, scripts, configurations, and test materials are described in the corresponding sections and folder READMEs.

```text
.
├── README.md
├── Documentation.pdf
├── LICENSE
├── .gitattributes
├── phenotype.tables/
│   ├── caas.discovery.groups.csv
│   ├── caas.results.csv
│   ├── human.equivalent.measurements.csv
│   ├── metadata.qualitative.traits.csv
│   ├── metadata.quantitative.traits.csv
│   ├── nhp.phenomic.dataset.csv
│   ├── trait.data.sources.csv
│   └── traits.redundancy.csv
├── genome.phenome.maps/
│   ├── ALL_RER_CONVERGE_UNIQUE_MAR23.txt
│   └── CAAS_FINAL_MAR23.tab
├── orthologous.alignments/
│   └── README.md
└── valenzuela2026.supplementary.information/
    ├── README.md
    ├── supplementary.tables/
    │   ├── Data_S1.xlsx
    │   ├── Data_S1.ods
    │   ├── Data_S2.xlsx
    │   ├── Data_S2.ods
    │   ├── Data_S3.xlsx
    │   └── Data_S3.ods
    ├── figure.creation/
    │   ├── MAIN/
    │   │   ├── InputFiles/
    │   │   └── Scripts/
    │   └── SUPP/
    │       ├── S1/ … S50/ (S11 is named s11/)
    │       └── S6,9,12,15,18,21_MERGED/
    └── genome.phenome.workflows/
        ├── caas/
        │   ├── inputs.and.settings/
        │   │   ├── alignments/
        │   │   ├── configs/
        │   │   └── total_background/
        │   └── fastcaas.fork/
        │       ├── README.md
        │       ├── fastcaas
        │       ├── Python modules and R permulation code
        │       ├── launch scripts
        │       ├── _tests/
        │       ├── _issues/
        │       └── _oldscript/
        └── rerconverge/
            └── link_tables_NOV21
```

## 1. Phenotypic data: `phenotype.tables/`

The [phenotype tables](phenotype.tables/) provide the phenotypic foundation of the archive. They connect species-level observations with trait definitions, quantitative summaries, categorical coverage, source references, and the species groups used for CAAS discovery. The CSV format supports direct inspection in spreadsheet applications and import into R, Python, or other analysis environments.

### Species-by-trait measurements

[`nhp.phenomic.dataset.csv`](phenotype.tables/nhp.phenomic.dataset.csv) contains the phenomic matrix. Its leading columns include `SpeciesBROAD` and `GroupName`, followed by trait identifiers such as `PTD00023` and `PTD00001`. This organization allows users to select species and traits for comparative analyses and connect the measurements to the accompanying metadata.

Use the trait identifiers to connect the matrix with the descriptive tables. Species labels provide the corresponding connection to trait-group definitions and sequence-based analyses. When combining tables, retain the supplied identifiers alongside descriptive trait names so that the provenance of each measurement remains traceable.

### Trait descriptions and sources

[`trait.data.sources.csv`](phenotype.tables/trait.data.sources.csv) records trait identifiers and names, descriptions, descriptions of values, source information, and source links. It is the starting point for understanding the biological meaning of an entry in the phenomic matrix and tracing the observations to the underlying references.

### Quantitative and qualitative metadata

[`metadata.quantitative.traits.csv`](phenotype.tables/metadata.quantitative.traits.csv) summarizes quantitative traits. Its fields include trait identifiers and names, species and family coverage, trait classes, and descriptive statistics: mean, standard deviation, minimum, quartiles, median, and maximum. These summaries help users assess the distribution and taxonomic coverage of a quantitative variable before selecting it for analysis.

[`metadata.qualitative.traits.csv`](phenotype.tables/metadata.qualitative.traits.csv) describes qualitative traits. It includes trait identifiers and names, covered species, trait types, category definitions, counts by category, and family coverage by category. It supports interpretation of categorical phenotypes and the composition of the categories compared in the study.

### Human-equivalent trait descriptions

[`human.equivalent.measurements.csv`](phenotype.tables/human.equivalent.measurements.csv) connects primate traits with human-equivalent measurements and terminology. The fields include measurement descriptions, source information, and mappings to human trait terms and identifiers from UMLS, EFO, causalDB, and the GWAS Catalog where recorded. These annotations provide a bridge between the comparative primate resource and human phenotype terminology.

### Discovery groups and trait redundancy

[`caas.discovery.groups.csv`](phenotype.tables/caas.discovery.groups.csv) describes the species assignments used to define CAAS discovery contrasts. Its columns include `ID`, `Species`, `TOP=1/BOTTOM=0`, and `Trait`. The group indicator records the top/bottom distinction used for the associated trait contrast.

[`traits.redundancy.csv`](phenotype.tables/traits.redundancy.csv) supplies the trait-redundancy curation table, with fields for identifiers, trait domains, trait names, and species and family counts. It provides context for interpreting related phenotypic measurements and the organization of the trait collection.

The folder also includes `caas.results.csv`. For the study's CAAS association tables and positional annotations, use the genome–phenome maps and the Data S2 workbook described below.

## 2. Molecular association results: `genome.phenome.maps/`

The [genome–phenome maps](genome.phenome.maps/) contain the principal result tables linking molecular evolutionary patterns to phenotypic traits. They represent two complementary approaches.

**Convergent amino acid substitutions (CAAS)** identify protein-alignment positions at which amino acid patterns are associated with the species groups defined by a phenotype. **Relative evolutionary rate (RER)** analyses investigate associations between gene evolutionary rates and phenotypic variation.

| File | Contents and role |
| --- | --- |
| [`CAAS_FINAL_MAR23.tab`](genome.phenome.maps/CAAS_FINAL_MAR23.tab) | CAAS results organized around genes, traits, alignment positions, amino acid patterns, and analysis scenarios. Use this table to explore the substitution-level component of the genome–phenome map. |
| [`ALL_RER_CONVERGE_UNIQUE_MAR23.txt`](genome.phenome.maps/ALL_RER_CONVERGE_UNIQUE_MAR23.txt) | RERconverge association results linking genes and traits, together with numerical analysis outputs. Use this table to explore the evolutionary-rate component of the genome–phenome map. |

The files are tab-delimited and retain their original analysis filenames. Data S2 in the supplementary materials provides related result tables and annotations, including CAAS contrast metadata, mutation results, mammalian-equivalent positions, and an RER results worksheet labeled with an FDR threshold of 0.05. Consult the workbook's contents index when choosing the result representation appropriate to a particular question.

## 3. Sequence alignments: `orthologous.alignments/`

The [alignment README](orthologous.alignments/README.md) describes the **PGA Primate 1-to-1 Ortholog Alignments**, distributed as the versioned dataset **`PGA-1to1-v1.0`**. The collection contains primate coding-sequence and protein alignments associated with Kuderna et al. (2023) and used in Valenzuela et al. (2026).

Version 1.0 has a dataset creation date of **28 October 2026**. The folder README provides the Google Drive download location, archive sizes, package organization, citation information, attribution, and dataset-specific license information.

The dataset is packaged in ZIP and TAR.GZ formats. Both archives contain the same collection; choose either format. The extracted package includes:

| Directory or file | Contents |
| --- | --- |
| `cds/` | 16,133 unfiltered coding-sequence alignments in FASTA format. |
| `cds.filtered/` | 16,133 filtered coding-sequence alignments in FASTA format. |
| `cds.filtered.translated/` | 16,133 protein alignments corresponding to the filtered coding sequences, in PHYLIP format. |
| `README.md` | Dataset description, filename conventions, references, and attribution. |
| `LICENSE` | Dataset-specific license notice and the full CC BY 4.0 legal text. |
| `md5.codes.txt` | Relative alignment paths and their MD5 checksums. |

The package contains **48,399 alignment files** and requires approximately **15.1 GB** after extraction. Sequence headers identify the taxa represented in each alignment. The species name included in an alignment filename should be interpreted together with those headers when identifying the sequences in a file.

## 4. Study materials: `valenzuela2026.supplementary.information/`

The [supplementary-materials README](valenzuela2026.supplementary.information/README.md) is the detailed code and data guide. It describes the supplementary datasets, figure-analysis organization, genome–phenome workflows, and file formats. It also provides an inventory with relative paths, MD5 checksums, and descriptions for the supplementary workbooks and figure-analysis files.

### 4.1. Supplementary workbooks: `supplementary.tables/`

The [supplementary tables](valenzuela2026.supplementary.information/supplementary.tables/) contain three datasets, each provided as an Excel workbook (`.xlsx`) and an OpenDocument Spreadsheet (`.ods`). The two formats are alternative representations of each dataset. Each workbook contains an index describing its worksheets.

| Dataset | Subject | Included information |
| --- | --- | --- |
| **Data S1** | Phenotypic data and trait organization | Non-human primate phenomic observations; trait classification and sources; human-equivalent measurements; quantitative and qualitative metadata; trait redundancy; allometric traits; and trait top/bottom selections. |
| **Data S2** | Genetic associations and positional annotations | Filtered gene-alignment metadata; CAAS contrast metadata; CAAS results; RER association results; CAAS mutation results; and mammalian-equivalent positions. |
| **Data S3** | Functional interpretation and trait-pair relationships | Permulation-based functional-enrichment tables; phylogenetic generalized least-squares (PGLS) enrichment results; and trait-pair overlap results. |

The standalone phenotype CSVs support direct access to individual tables, while these workbooks collect the study's supplementary datasets into structured, multi-worksheet documents.

### 4.2. Figure analyses: `figure.creation/`

The [figure-creation directory](valenzuela2026.supplementary.information/figure.creation/) contains R scripts and the input materials used in the study's figure analyses. The organization follows the distinction between main-text figures and supplementary figures.

- **`MAIN/InputFiles/`** holds shared inputs for main-figure analyses, including association results, trait descriptions, trait-group configurations, phylogenetic information, annotations, serialized R objects, and plotting palettes.
- **`MAIN/Scripts/`** contains 13 R scripts associated with main figures and individual figure panels. Their filenames identify the relevant analyses and outputs.
- **`SUPP/S1/` through `SUPP/S50/`** contain the individual supplementary-figure analyses. The folder for Figure S11 is named `s11/`. The collection includes 50 individual supplementary-figure scripts and figure-specific auxiliary files.
- **`SUPP/S6,9,12,15,18,21_MERGED/`** contains six additional R scripts for the corresponding merged analyses.

Together, these directories contain **69 R scripts**. The scripts document data loading, filtering, joins, statistical calculations, and plotting choices. Auxiliary files include text tables, trait configurations, phylogenetic trees, chromosome-specific conservation and variation information, and serialized R objects. Some analyses save R objects with `saveRDS()`; the reading and saving commands in each script specify how those files are used.

Use the detailed supplementary README to locate inputs and identify the scripts associated with a particular figure. Its checksum inventory covers the six workbook files and 350 files in the figure-analysis collection.

### 4.3. Genome–phenome analyses: `genome.phenome.workflows/`

The [genome–phenome workflows](valenzuela2026.supplementary.information/genome.phenome.workflows/) preserve the code, inputs, and settings associated with the study's CAAS and RERconverge analyses. This collection documents the computational context behind the association results and complements the figure-analysis scripts.

#### CAAS inputs and settings: `caas/inputs.and.settings/`

This directory groups the materials used to define and run the phenotype-specific CAAS analyses:

- **`alignments/`** contains protein alignments supplied as inputs to CAAS discovery.
- **`configs/`** contains phenotype-specific configurations and permulation values, organized into qualitative and quantitative trait categories. These files record the trait and species-group information used by the analysis workflows.
- **`total_background/`** contains `BG_QUALITATIVE.txt` and `BG_QUANTITATIVE.txt`, which list trait–gene background associations for the corresponding analyses.

#### Archived implementation: `caas/fastcaas.fork/`

The archived implementation includes the `fastcaas` executable, Python analysis modules, R permulation code, example launch scripts, and its own README. The `_tests/` directory holds test and example materials; `_issues/` and `_oldscript/` retain development materials and earlier scripts associated with the archived implementation.

**CAAStools is called `fastcaas` in this archive because it is a version optimized to process many phenotypes within a single workflow. It is included for reproducibility and to document the implementation used in the study. We advise against running this archived version. For further analyses, use [CAAStools](https://github.com/linudz/caastools) and follow its documentation.**

Additional scripts for running the analyses on a computing cluster are available on request from [fabio.barteri@upf.edu](mailto:fabio.barteri@upf.edu).

#### RERconverge materials: `rerconverge/`

The file [`link_tables_NOV21`](valenzuela2026.supplementary.information/genome.phenome.workflows/rerconverge/link_tables_NOV21) records the Google Drive location for the RERconverge materials. Follow that location to access the externally distributed collection associated with the evolutionary-rate analyses.

## 5. Accessing and working with the archive

### Browsing the data

For phenotype definitions and observations, start with `phenotype.tables/`. For molecular associations, use `genome.phenome.maps/` together with Data S2. For functional-enrichment and trait-overlap analyses, consult Data S3. For the sequence collection, follow the versioned download instructions in `orthologous.alignments/README.md`. For figure-specific analyses and their inputs, use the supplementary-materials README and `figure.creation/`.

### Downloading Git LFS files

Large files tracked through Git Large File Storage are identified in [`.gitattributes`](.gitattributes). Use a Git LFS-enabled checkout to retrieve these data:

```bash
git lfs install
git clone https://github.com/pgarchive/data.git
cd data
git lfs pull
```

The two principal genome–phenome maps and selected shared figure inputs use Git LFS. The full versioned alignment package is accessed through the download instructions in its dedicated README.

### Connecting the different materials

Trait identifiers and trait names connect phenotypic observations, metadata, discovery groups, and association results. Gene identifiers connect molecular associations to alignments and functional annotations. Alignment positions and the positional annotations in Data S2 support interpretation of individual CAAS results.

Retain the original identifiers, filenames, and dataset version when extracting subsets or documenting an analysis. Use the source-reference table for phenotypic provenance, the workbook indexes for worksheet organization, and the relevant configuration files for phenotype-specific contrast definitions.

### Documentation and integrity information

[`Documentation.pdf`](Documentation.pdf) provides additional documentation of the PGA data resource. The [supplementary README](valenzuela2026.supplementary.information/README.md) provides the detailed directory guide and MD5 inventory for supplementary tables and figure-analysis files. The alignment package supplies its own `md5.codes.txt` file, with paths relative to the extracted dataset directory.

## 6. Licenses and attribution

The repository includes the [GNU General Public License v3.0](LICENSE). Consult that license and the notices supplied with individual components when reusing repository materials.

The separately distributed **`PGA-1to1-v1.0` alignment dataset** has its own **Creative Commons Attribution 4.0 International (CC BY 4.0)** license. Its owner and licensor is **Universitat Pompeu Fabra (UPF)**. The dataset-specific notice and full legal text are included in the alignment package's `LICENSE` file and described in the [alignment README](orthologous.alignments/README.md).

The resource is associated with the **Institute of Evolutionary Biology (IBE, CSIC–UPF), Barcelona, Spain**. When using the data, cite the relevant study and retain the attribution and license information accompanying the materials used.

## 7. References and contact

- Valenzuela et al. (2026). *A phylogenetic protein-coding genome-phenome map of complex traits across 224 primate species*. **bioRxiv**, version 3. [doi:10.1101/2025.09.08.674744](https://doi.org/10.1101/2025.09.08.674744).
- Kuderna et al. (2023). *A global catalog of whole-genome diversity from 233 primate species*. **Science**, 380, 906–913. [doi:10.1126/science.abn7829](https://doi.org/10.1126/science.abn7829).

For questions about the archived workflows or requests for additional cluster execution scripts, contact **Dr. Fabio Barteri** at [fabio.barteri@upf.edu](mailto:fabio.barteri@upf.edu).
