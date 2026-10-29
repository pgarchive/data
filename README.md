# The Primate Genome-phenome Archive (PGA)

This repository contains core data files for the **Primate Genome-phenome Archive (PGA)**, a curated resource integrating genetic and phenotypic information from non-human primates.

The PGA dataset is part of the project *Valenzuela et al., 2025*, and provides a foundation for evolutionary, biomedical, and comparative analyses linking genotype to phenotype across primates.

## Repository Structure
```
├── genetic.information
│   ├── ALL_RER_CONVERGE_UNIQUE_MAR23.txt
│   └── CAAS_FINAL_MAR23.tab
├── LICENSE
├── phenotype.tables
│   ├── caas.discovery.groups.csv
│   ├── caas.results.csv
│   ├── human.equivalent.measurements.csv
│   ├── metadata.qualitative.traits.csv
│   ├── metadata.quantitative.traits.csv
│   ├── nhp.phenomic.dataset.csv
│   ├── trait.data.sources.csv
│   └── traits.redundancy.csv
└── README.md
```

### `genetic.information/`

This folder contains gene-level evolutionary signals, including:

- **CAAS_FINAL_MAR23.tab**: Final list of convergent amino acid substitutions (CAAS) detected across primates.
- **ALL_RER_CONVERGE_UNIQUE_MAR23.txt**: Genes showing significant rate shifts (RERconverge) associated with phenotype convergence.

### `phenotype.tables/`

This folder includes phenotype tables describing ~230 traits measured across ~220 primate species. The data include:

- Trait metadata (quantitative and qualitative)
- Cross-species harmonization (e.g., human equivalent units)
- Trait grouping and redundancy curation
- Results of gene-trait association analyses

For a detailed explanation of each file, please refer to the official PGA documentation:

📄 [Download the PDF documentation](https://github.com/pgarchive/data/blob/main/Documentation.pdf)  

## License

This repository is distributed under the **GNU General Public License v3.0**.  
See the [LICENSE](LICENSE) file for details.

---

*Part of the project "Valenzuela et al., 2025".*

## Supplementary materials

The [supplementary materials](valenzuela2026.supplementary.information/) contain the supplementary datasets and figure-analysis scripts. Details on file organization, available datasets, execution limitations, and the complete file inventory with MD5 checksums are provided in the code and data guide ([Markdown](valenzuela2026.supplementary.information/code.and.data.guide.md) | [PDF](valenzuela2026.supplementary.information/code.and.data.guide.pdf)).
