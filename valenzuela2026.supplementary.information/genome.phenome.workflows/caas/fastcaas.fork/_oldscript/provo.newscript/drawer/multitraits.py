#!/usr/bin/env python3

# Fastcaas - fc_lim_indexed_traits.py module.
# Scan alignments with multiple traits.

from index_configs import *
from slice_and_index import *


### SETTINGS

alignments = "alignments"
the_traits_folder = "lesscfgs"

### FUNCTIONS

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
    
    return z


# FUNCTION check_scenario()
# checks the scenario

def check_scenario(input_string):
    scenario = "4"
    twosides = input_string.split("/")
    fg = len(list(twosides[0]))
    bg = len(list(twosides[0]))

    if fg == 1 and bg == 1:
        scenario = "1"
    
    elif fg == 1:
        scenario = "2"
    elif bg == 1:
        scenario = "3"
    
    return scenario

# FUNCTION fetch_caas():
# fetches caas per each thing

def fetch_caas(genename, processed_position, list_of_traits, output_file, admitted_scenarios = ["1","2","3"], maxgaps_fg = 0, maxgaps_bg = 0):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    positive_traits = []

    # Test 1: is this position a caas?:
    
    aas = list(processed_position.aas2traits.keys())

    for a in aas:
    
        for trait in valid_traits:
            fg = trait + "_1"
            bg = trait + "_0"
    
            if a != "-":

                if fg in processed_position.aas2traits[a] or bg in processed_position.aas2traits[a]:
                    positive_traits.append(trait)

                if fg in processed_position.aas2traits[a] and bg in processed_position.aas2traits[a]:
                    positive_traits.remove(trait)
            
            else:
                if a.count(fg) > maxgaps_fg:
                    if trait in positive_traits:
                        positive_traits.remove(trait)
                if a.count(bg) > maxgaps_bg:
                    if trait in positive_traits:
                        positive_traits.remove(trait)
    

    positive_traits = list(set(positive_traits))
    output_traits = []
    for x in positive_traits:
        aa_tag_fg_list = processed_position.trait2aas_fg[x]
        aa_tag_fg_list.sort()
        aa_tag_fg = "".join(aa_tag_fg_list)

        aa_tag_bg_list = processed_position.trait2aas_bg[x]
        aa_tag_bg_list.sort()
        aa_tag_bg = "".join(aa_tag_bg_list)

        tag = "/".join([aa_tag_fg, aa_tag_bg])

        scenario = check_scenario(tag)

        if scenario in admitted_scenarios:
            output_traits.append(x + "@" + tag)
            print("yes")

    if len(output_traits) > 0:
        out = open(output_file, "a")
        print(genename, processed_position.position, ",".join(output_traits), file=out)
        print(genename, processed_position.position, ",".join(output_traits))
        out.close()

### test

def scan_alignment_multi(alignment, the_traits_folder):
    p = slice(alignment)
    mcf = load_cfg_dictionary(the_traits_folder)
    for position in p.d:
        pro_position = process_position(position, mcf)
        fetch_caas(p.genename, pro_position, mcf.alltraits, "output.txt")


for a in glob.glob(alignments + "/*")[:10]:
    print(a)
    scan_alignment_multi(a, the_traits_folder)