### [The Primate Genome-Phenome Archive](../) / orthologous.alignments

# One to One Ortholog CDS Alignments

**Version: 1.0**

**Creation date: 28 October 2026**

This collection contains nucleotide and protein alignments of primate one-to-one orthologs associated with the work of Kuderna et al. (2023) and used in Valenzuela et al. (2026).

## Access version 1.0

**Dataset identifier:** `PGA-1to1-v1.0`  
**Download location (Google Drive):** [https://drive.google.com/drive/folders/1mprhR0bt___WG0Ue0Lh7B7Fyxf6qPZXF?usp=sharing](https://drive.google.com/drive/folders/1mprhR0bt___WG0Ue0Lh7B7Fyxf6qPZXF?usp=sharing)


Version 1.0 is packaged in two alternative formats:

| Archive | Size | Contents |
| --- | ---: | --- |
| `PGA-1to1-v1.0.zip` | 583,979,209 bytes (approximately 584 MB) | Complete version 1.0 dataset |
| `PGA-1to1-v1.0.tar.gz` | 580,191,755 bytes (approximately 580 MB) | Complete version 1.0 dataset |

Once access is available, open the download folder, select version 1.0, and download **one** of these archives. Both formats contain the same dataset; downloading both is unnecessary. Extract the archive to obtain the `PGA-1to1-v1.0/` directory. Allow approximately 15.1 GB of space for the extracted files, in addition to the downloaded archive.

## Archive contents

Each archive contains the following directory:

```text
PGA-1to1-v1.0/
├── README.md
├── LICENSE
├── md5.codes.txt
├── cds/
├── cds.filtered/
└── cds.filtered.translated/
```

| Directory | Contents | Format | Alignment files |
| --- | --- | --- | ---: |
| `cds/` | Unfiltered coding sequence (CDS) alignments | FASTA (`.fasta`) | 16,133 |
| `cds.filtered/` | Filtered CDS alignments | FASTA (`.fasta`) | 16,133 |
| `cds.filtered.translated/` | Protein alignments corresponding to the filtered CDS | PHYLIP (`.phy`) | 16,133 |

There are **48,399 alignment files** in total. The included `README.md` describes the collections, filename conventions, references, and attribution. The `md5.codes.txt` file provides the MD5 checksums of all alignment files as a tab-separated table with `path` and `md5` columns; paths are relative to the extracted dataset directory.

The species name in a filename does not describe the full taxonomic composition of an alignment. Sequence headers identify the taxa represented within each file.

## License and attribution

**Owner and licensor:** Universitat Pompeu Fabra (UPF)  
**Copyright:** © 2026 Universitat Pompeu Fabra (UPF).

The dataset **PGA-1to1-v1.0** is distributed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). Copying, redistribution, and adaptation, including commercial use, are permitted under its terms. When sharing the dataset or modified versions, give appropriate credit, retain the supplied notices, link to the license, and indicate any modifications.

The dataset-specific license notice and the full legal text are included in the **`LICENSE` file inside each archive**.

**Attribution:** Universitat Pompeu Fabra (UPF). *PGA-1to1-v1.0* (2026). The Primate Genome-Phenome Archive.

Please also cite the original publications listed below.

## Version history

| Version | Dataset creation date | Dataset identifier | Changes | Access status |
| --- | --- | --- | --- | --- |
| 1.0 | 28 October 2026 | `PGA-1to1-v1.0` | Initial packaged version: unfiltered CDS, filtered CDS, and corresponding protein alignments; README, CC BY 4.0 license, and alignment MD5 checksums included. | Google Drive folder not yet public. |

This history will be updated for subsequent versions with the release date, dataset identifier, changes, and download location. Version identifiers distinguish the distributed packages; the date above is the dataset creation date, not a public release date.

## References

- Kuderna et al. (2023). *A global catalog of whole-genome diversity from 233 primate species*. **Science**, 380, 906–913. [doi:10.1126/science.abn7829](https://doi.org/10.1126/science.abn7829).
- Valenzuela et al. (2026). *A phylogenetic protein-coding genome-phenome map of complex traits across 224 primate species*. **bioRxiv**, version 3. [doi:10.1101/2025.09.08.674744](https://www.biorxiv.org/content/10.1101/2025.09.08.674744v3).

## Dataset attribution and contact

This dataset is attributed to the [Institute of Evolutionary Biology (IBE, CSIC–UPF)](https://www.ibe.upf-csic.es/), Barcelona, Spain.

**Contact:** Dr. Fabio Barteri — [fabio.barteri@upf.edu](mailto:fabio.barteri@upf.edu).
