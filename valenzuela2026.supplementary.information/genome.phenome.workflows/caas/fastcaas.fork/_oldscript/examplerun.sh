#!/bin/bash

# Sample run for FAST CAAS. Just run
#   bash samplerun.sh       or
#   ./samplerun.sh

# You might need to setup your python environment first

# On fasta files
./fastcaas.py -c examples/examplebinary.cfg -f examples/80alis -o examples/exampleresult.tab

# On phylip files
./fastcaas.py -c examples/examplebinary.cfg -f examples/80alis.phylip --fmt phylip -o examples/exampleresult.phylip.tab

# On phylip-relaxed files
./fastcaas.py -c examples/examplebinary.cfg -f examples/80alis.phylip-relaxed --fmt phylip-relaxed -o examples/exampleresult.phylip-relaxed.tab


# On clustal files
./fastcaas.py -c examples/examplebinary.cfg -f examples/80alis.clustal --fmt clustal -o examples/exampleresult.clustal.tab

# On stockholm files
./fastcaas.py -c examples/examplebinary.cfg -f examples/80alis.stockholm --fmt stockholm -o examples/exampleresult.stockholm.tab
