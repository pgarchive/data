#!/bin/bash +x

/bin/date

# For each gene, launch fastcaas. If the -j option is the same, the same folder will be updated with new results

# WARNING. This is the entrance to the rabbit hole

# Generic usage:
#   for alignment in d:
#       ./fastcaas_alt -a $alignment -t $traitfile --species #speciesfile -j $jobname(the same for all files)


# Init mode. Run the initialisation.

./fastcaas init -j _issues/rabbithole/_simulate -t _tests/BodyMass_kg_permulation.cfg --species _tests/_samplerun/sp2fam.210727.tab --bootstrap random --phylogeny _issues/permulation_implementation/inputs/traits.speciestree.210512.nh --traitvals _tests/BodyMass_kg_permulation.txt

# Run mode. Example on ZACN file...
./fastcaas run -a _tests/ZACN.Homo_sapiens.filter2.fasta -t _tests/BodyMass_kg_permulation.cfg --species _tests/_samplerun/sp2fam.210727.tab -j _issues/rabbithole/_simulate

# Finally, run the collate mode to merge the results into the bootstrap

./fastcaas collate -j _issues/rabbithole/_simulate

/bin/date