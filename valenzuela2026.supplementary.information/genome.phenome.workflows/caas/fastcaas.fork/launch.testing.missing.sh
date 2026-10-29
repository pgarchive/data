#!/bin/bash +x

/bin/date

# For each gene, launch fastcaas. If the -j option is the same, the same folder will be updated with new results

# Generic usage:
#   for alignment in d:
#       ./fastcaas_alt -a $alignment -t $traitfile --species #speciesfile -j $jobname(the same for all files)


# Init mode. Run the initialisation.

./fastcaas init -j _tests/_samplerun/_results -t _tests/BodyMass_kg_permulation.cfg.wspecies --species _tests/_samplerun/sp2fam.210727.tab.wspecies --bootstrap permulations --phylogeny _issues/permulation_implementation/inputs/traits.speciestree.210512.nh --traitvals _tests/BodyMass_kg_permulation.txt

# Run mode. Example on two files...
./fastcaas run -a _tests/_samplerun/A1BG.Homo_sapiens.filter2.phy -t _tests/BodyMass_kg_permulation.cfg.wspecies --species _tests/_samplerun/sp2fam.210727.tab.wspecies -j _tests/_samplerun/_results --max_miss 3 --max_gaps 0
./fastcaas run -a _tests/_samplerun/ZSWIM8.Homo_sapiens.filter2.fasta -t _tests/BodyMass_kg_permulation.cfg.wspecies --species _tests/_samplerun/sp2fam.210727.tab.wspecies -j _tests/_samplerun/_results --max_miss 3 --max_gaps 0

# Finally, run the collate mode to merge the results into the bootstrap

./fastcaas collate -j _tests/_samplerun/_results

/bin/date
