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

MODULE NAME: dphase.py
DESCRIPTION: runs the caas on one single alignment
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: fascaas pipeline script

'''

from init_bootstrap import *
from disco import process_position
from caas_id import filter_for_gaps, iscaas
from alimport import *
from os.path import exists



# FUNCTION fetch_caas():
# fetches caas per each thing

def caasboot(genename, processed_position, list_of_traits, output_file, maxgaps_fg, maxgaps_bg, maxgaps_all, admitted_scenarios = ["1","2","3"], cycles = 1000):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    # Filter for the number of gaps

    output_traits = []

    for trait in valid_traits:
        if filter_for_gaps(max_bg = maxgaps_bg, max_fg = maxgaps_fg, max_all = maxgaps_all, gfg = processed_position.trait2gaps_fg[trait], gbg = processed_position.trait2gaps_bg[trait]) == True:
            # Filter for scenario

            aa_tag_fg_list = processed_position.trait2aas_fg[trait]
            aa_tag_fg_list.sort()
            aa_tag_fg = "".join(aa_tag_fg_list)

            aa_tag_bg_list = processed_position.trait2aas_bg[trait]
            aa_tag_bg_list.sort()
            aa_tag_bg = "".join(aa_tag_bg_list)

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag)

            if check.caas == True and check.scenario in admitted_scenarios:
                output_traits.append(trait)

    # Print the output

    if exists(output_file):
        out = open(output_file, "a")
    else:
        out = open(output_file, "w")
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))
    pvalue = str(len(output_traits)/cycles)
    traitline = ",".join(output_traits)

    outline = "\t".join([position_name, count, pvalue, traitline])

    print(outline, file=out)
        


    out.close()

# FUNCTION disco_bootstrap()
# Launches the bootstrap in several lines

def disco_bootstrap(input_alignment, input_cfg, input_tree, max_fg_gaps, max_bg_gaps, max_overall_gaps, admitted_scenarios, output_file):

    # Step 0: global_variables

    genename = input_alignment.split("/")[1].split(".")[0]

    # Step 1: import the trait into a trait object (load_cfg from pindex.py)
    trait_object = simtrait(input_cfg, input_tree)

    # Step 2: import the alignment int a processed position object (slice from alimport.py)
    p = slice(input_alignment)

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = trait_object), p.d)

    # Step 4: extract the raw caas
    for position in processed_positions:
        caasboot( genename,
                    position,
                    trait_object.alltraits,
                    maxgaps_bg= max_bg_gaps,
                    maxgaps_fg= max_fg_gaps,
                    maxgaps_all= max_overall_gaps,
                    admitted_scenarios=admitted_scenarios,
                    output_file = output_file
                    )


### Function assign_pvalue

def assign_pvalue(bootstrapfile):
    return 0

# Test
'''
tree_path = "bt_testing/traits.speciestree.210512.nh"
cfg_path = "tests/AOB_config.txt"
alignment = "tests/ZACN.Homo_sapiens.filter2.fasta"

imported_bootstrap = simtrait(cfg_path, tree_path)

disco_bootstrap(    input_alignment = alignment,
                    input_cfg = cfg_path,
                    input_tree = tree_path,
                    max_fg_gaps = 0,
                    max_bg_gaps = 0,
                    max_overall_gaps = 0,
                    admitted_scenarios = ["1","2","3"],
                    output_file = "provaboot.tuf")
'''