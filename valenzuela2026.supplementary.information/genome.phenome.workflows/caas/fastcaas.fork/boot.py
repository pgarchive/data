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
from caas_id import filter_for_gaps, filter_for_missings, iscaas
from alimport import *

from os.path import exists
import functools


# FUNCTION fetch_caas():
# fetches caas per each thing



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, admitted_scenarios = ["1","2","3"], cycles = 1000):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    def filter_trait(trait, the_processed_position, max_bg_gaps, max_fg_gaps, max_all_gaps, max_bg_miss, max_fg_miss, max_all_miss):
        
        the_gfg = the_processed_position.trait2gaps_fg[trait]
        the_gbg = processed_position.trait2gaps_bg[trait]

        the_mfg = the_processed_position.trait2miss_fg[trait]
        the_mbg = processed_position.trait2miss_bg[trait]


        if filter_for_gaps(max_bg = max_bg_gaps, max_fg = max_fg_gaps, max_all = max_all_gaps, gfg = the_gfg, gbg = the_gbg) == False:
            return False

        elif filter_for_missings(max_m_bg = max_bg_miss, max_m_fg = max_fg_miss, max_m_all = max_all_miss, mfg = the_mfg, mbg = the_mbg) == False:
            return False

        else:
            aa_tag_fg_list = the_processed_position.trait2aas_fg[trait]
            aa_tag_fg_list.sort()
            aa_tag_fg = "".join(aa_tag_fg_list)

            aa_tag_bg_list = the_processed_position.trait2aas_bg[trait]
            aa_tag_bg_list.sort()
            aa_tag_bg = "".join(aa_tag_bg_list)

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag)

            if check.caas == True and check.scenario in admitted_scenarios:
                return True
    
    output_traits = filter(functools.partial(
        filter_trait,
        the_processed_position = processed_position,
        max_bg_gaps = maxgaps_bg,
        max_fg_gaps = maxgaps_fg,
        max_all_gaps = maxgaps_all,

        max_bg_miss = maxmiss_bg,
        max_fg_miss = maxmiss_fg,
        max_all_miss = maxmiss_all,
    ),
                    valid_traits)



    output_traits = list(output_traits)

    # Return the line
    
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(output_traits))
    pvalue = str(len(output_traits)/cycles)
    traitline = ",".join(output_traits)

    outline = "\t".join([position_name, count, pvalue, traitline])

    return outline
        

# FUNCTION disco_bootstrap()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, simulated_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_scenarios, output_file):


    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = simulated_traits, species_in_alignment = sliced_object.species), sliced_object.d)
    the_genename = sliced_object.genename

    # Step 4: extract the raw caas

    output_lines = map(
        functools.partial(
            caasboot,
            list_of_traits = simulated_traits.alltraits,
            genename = the_genename,
            maxgaps_fg = max_fg_gaps,
            maxgaps_bg = max_bg_gaps,
            maxgaps_all = max_overall_gaps,

            maxmiss_fg = max_fg_miss,
            maxmiss_bg = max_bg_miss,
            maxmiss_all = max_overall_miss,

            admitted_scenarios = the_admitted_scenarios,
            cycles = 1000) ,processed_positions
    )

    output_lines = list(output_lines)

    ooout = open(output_file, "w")

    for line in output_lines:
        print(line + "\t" + trait_config_file, file=ooout)
    
    ooout.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d

# Test

'''
species_path = "tests/sp2fam.210727.tab"
cfg_path = "tests/AOB_config.txt"
alignment = "tests/ZACN.Homo_sapiens.filter2.fasta"

print("Executing random bootstrap")
imported_bootstrap = simtrait(cfg_path, species_path, mode = "random")
imported_bootstrap.print_traits("provaboot.ZACN.random.traitlist.tuf")

disco_bootstrap(    input_alignment = alignment,
                    simulated_traits = imported_bootstrap,
                    max_fg_gaps = 0,
                    max_bg_gaps = 0,
                    max_overall_gaps = 0,
                    admitted_scenarios = ["1","2","3"],
                    output_file = "provaboot.ZACN.random.result.tuf")


print("Executing phylogeny-guided bootstrap")
imported_bootstrap = simtrait(cfg_path, species_path, mode = "phylogeny-guided")
imported_bootstrap.print_traits("provaboot.ZACN.phy.traitlist.tuf")

disco_bootstrap(    input_alignment = alignment,
                    simulated_traits = imported_bootstrap,
                    max_fg_gaps = 0,
                    max_bg_gaps = 0,
                    max_overall_gaps = 0,
                    admitted_scenarios = ["1","2","3"],
                    output_file = "provaboot.ZACN.phy.result.tuf")
'''