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

MODULE NAME:  caas_id.py
DESCRIPTION:  CAAS identification form imported alignment and indexed phenotypes
DEPENDENCIES: alimport.py and pindex.py
CALLED BY:    fastcaas_core.py

TABLE OF CONTENTS
------------------------------------------

process_position()          processes a position from an imported alignment.
                            The output will be a dictionary that points each 
                            aminoacid (gaps included) to the
                            species sharing it.

check_scenario()            checks the scenario

fetch_caas()                fetches caas per each thing
'''                                                       

from pindex import *
from alimport import *
from os.path import exists

# Function process_position()
# processes a position from an imported alignment. The output will be
# a dictionary that points each aminoacid (gaps included) to the
# species sharing it.

def process_position(position, multiconfig):

    class processed_position():
        def __init__(self):
            self.position = ""
            self.aas2species = {}
            self.aas2traits = {}
            self.trait2aas_fg = {}
            self.trait2aas_bg = {}

            self.trait2ungapped_fg = {}
            self.trait2ungapped_bg = {}

            self.trait2gaps_fg = {}
            self.trait2gaps_bg = {}
            self.trait2gaps_all = {}

            self.gapped = []
        
    z = processed_position()

    # Load aas2species
    for x in position.keys():
        z.position = position[x].split("@")[1]
        try:
            z.aas2species[position[x].split("@")[0]].append(x)
        except:
            z.aas2species[position[x].split("@")[0]] = [x]

    # Load aas2traits    
    for key in z.aas2species.keys():
        traits = []
        for species in z.aas2species[key]:
            try:
                for v in multiconfig.s2t[species]:
                    if v not in traits:
                        traits.append(v)
            except:
                pass
        
        z.aas2traits[key] = traits

        for t in traits:
            if t[-2:] == "_1" and key != "-":
                try:
                    z.trait2aas_fg[t[:-2]].append(key)
                except:
                    z.trait2aas_fg[t[:-2]] = [key]
                    pass
            if t[-2:] == "_0" and key != "-":
                try:
                    z.trait2aas_bg[t[:-2]].append(key)
                except:
                    z.trait2aas_bg[t[:-2]] = [key]
                    pass
    
    try:
        z.gapped = z.aas2species["-"]
    except:
        pass

    # Determine Ungapped Species

    for trait in z.trait2aas_bg.keys():
        
        # Ungapped species

        ugfg = list(set(multiconfig.trait2fg[trait]) - set(z.gapped))
        z.trait2ungapped_fg[trait] = ugfg

        ugbg = list(set(multiconfig.trait2bg[trait]) - set(z.gapped))
        z.trait2ungapped_bg[trait] = ugbg

        # Number of gaps

        gfg = len(set(multiconfig.trait2fg[trait]).intersection(set(z.gapped)))
        gbg = len(set(multiconfig.trait2bg[trait]).intersection(set(z.gapped)))
        gall = gfg + gbg

        z.trait2gaps_fg[trait] = gfg
        z.trait2gaps_bg[trait] = gbg
        z.trait2gaps_all[trait] = gall
    return z


# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True
    all_g = gfg + gbg
    
    if all_g > max_all:
        out = False

    if gfg > max_fg:
        out = False

    if gbg > max_bg:
        out = False

    return out

# FUNCTION check_scenario()
# checks the scenario

def iscaas(input_string):
    
    class caaspositive():
        def __init__(self):
            self.caas = True
            self.scenario = "4"
    
    z = caaspositive()

    twosides = input_string.split("/")
    fg = list(twosides[0])
    bg = list(twosides[1])

    # Is this a CAAS?
    for x in fg:
        if x in bg:
            z.caas = False
            break

    # What is the scenario?
    if len(fg) == 1 and len(bg) == 1:
        z.scenario = "1"
    
    elif len(fg) == 1:
        z.scenario = "2"
    elif len(bg) == 1:
        z.scenario = "3"
    
    if len(fg) == 0 or len(bg) == 0:
        z.scenario = "null"
    
    return z

# FUNCTION fetch_caas():
# fetches caas per each thing

def fetch_caas(genename, processed_position, list_of_traits, output_file, maxgaps_fg, maxgaps_bg, maxgaps_all, admitted_scenarios = ["1","2","3"]):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    # Filter for the number of gaps

    for trait in valid_traits:
        if filter_for_gaps(max_bg = maxgaps_bg, max_fg = maxgaps_fg, max_all = maxgaps_all, gfg = processed_position.trait2gaps_fg[trait], gbg = processed_position.trait2gaps_bg[trait]) == False:
                valid_traits.remove(trait)

    output_traits = []

    # Filter for scenario
    for x in valid_traits:
        aa_tag_fg_list = processed_position.trait2aas_fg[x]
        aa_tag_fg_list.sort()
        aa_tag_fg = "".join(aa_tag_fg_list)

        aa_tag_bg_list = processed_position.trait2aas_bg[x]
        aa_tag_bg_list.sort()
        aa_tag_bg = "".join(aa_tag_bg_list)

        tag = "/".join([aa_tag_fg, aa_tag_bg])

        check = iscaas(tag)

        if check.caas == True and check.scenario in admitted_scenarios:
            output_traits.append(x + "@" + tag + "@scenario" + check.scenario)
    

    # Print the output

    if len(output_traits) > 0:

        if exists(output_file):
            out = open(output_file, "a")
        else:
            out = open(output_file, "w")

        for trait in output_traits:

            traitname = trait.split("@")[0]            
            change = trait.split("@")[1]
            thescenario = trait.split("@")[2]

            fg_species_number = str(len(processed_position.trait2ungapped_fg[traitname]))
            bg_species_number = str(len(processed_position.trait2ungapped_bg[traitname]))

            print(  "\t".join(
                [genename,
                    traitname,
                    processed_position.position,
                    change,
                    thescenario,
                    fg_species_number,
                    bg_species_number,
                    str(processed_position.trait2gaps_fg[traitname]),
                    str(processed_position.trait2gaps_bg[traitname]),
                    ",".join(processed_position.trait2ungapped_fg[traitname]),
                    ",".join(processed_position.trait2ungapped_bg[traitname])]
            ), file = out)

        out.close()