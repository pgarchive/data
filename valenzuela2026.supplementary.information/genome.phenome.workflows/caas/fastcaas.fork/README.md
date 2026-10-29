# FASTCAAS v 0.7

## Background
The identification of the molecular mechanisms that subtend convergent phenotypic adaptation is a prime interest for researchers working in evolutionary genomics. In a set of orthologous proteins, Convergent Amino Acid Substitutions (CAAS) will occur in those positions in which an amino acid converges among species that share a similar phenotype. They can emerge as an effect of positive selection and represent an interesting lead to unravel the molecular mechanisms underlying phenotypic variation.

**FASTCAAS** is a python tool that identifies CAAS from Multiple Sequence Alignments (MSA) of orthologous proteins. Given two groups of species, fastcaas will process one or more MSAs to identify all those positions in which the amino acids differ between the two groups. 

**FASTCAAS** is a command-line tool written in Python 3 and relies on the biopython library. It is compatible with different platforms and provided with parallelization options to process a high number of MSA input files.


**WARNING** We're still working on updated documentation.