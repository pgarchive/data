'''                                                          
,------. ,---.   ,---. ,--------. ,-----.  ,---.    ,---.   ,---.   
|  .---'/  O  \ '   .-''--.  .--''  .--./ /  O  \  /  O  \ '   .-'  
|  `--,|  .-.  |`.  `-.   |  |   |  |    |  .-.  ||  .-.  |`.  `-.  
|  |`  |  | |  |.-'    |  |  |   '  '--'\|  | |  ||  | |  |.-'    | 
`--'   `--' `--'`-----'   `--'    `-----'`--' `--'`--' `--'`-----'                                                                  

FAST Convergent Amino Acid Substitution identification tool

Authors:    Fabio Barteri (fabio.barteri@upf.edu),
            Xavier Farré (xfarrer@igtp.cat),
            David de Juan (david.juan@upf.edu).

MODULE NAME: disco.py
DESCRIPTION: runs the caas discovery on one single alignment. Returns non-validated caas candidate positions.
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: pipeline manager.

'''


from caas_id import *
from alimport import *
from pindex import *
import os

### FUNCTION discovery()
### Scans one single alignment to identify the CAAS

def discovery(input_cfg, sliced_object, pvalue_dictionary, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, admitted_scenarios, output_file):

    # Step 1: import the trait into a trait object (load_cfg from pindex.py)
    trait_object = load_cfg(input_cfg)

    # Step 2: import the alignment int a processed position object (slice from alimport.py)
    p = sliced_object

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = trait_object, species_in_alignment = p.species), p.d)

    # Step 4: extract the raw caas
    for position in processed_positions:
        fetch_caas( p.genename,
                    position,
                    trait_object.alltraits,
                    pvdict= pvalue_dictionary,

                    maxgaps_bg= max_bg_gaps,
                    maxgaps_fg= max_fg_gaps,
                    maxgaps_all= max_overall_gaps,

                    maxmiss_bg= max_bg_miss,
                    maxmiss_fg= max_fg_miss,
                    maxmiss_all= max_overall_miss,

                    admitted_scenarios=admitted_scenarios,
                    output_file = output_file
                    )