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

MODULE NAME: INIT BOOTSTRAP
DESCRIPTION: Initialises the bootstrap through different strategy (random, phylogeny and permulations)
DEPENDENCIES: none
CALLED BY: boot.py
'''

import dendropy, random

# FUNCTION simtrait() Simulate trait function

def simtrait(trait_config_file, tree_file, output_prefix  ="none", mode = "random", cycles = 1000, tree_format = "newick"):
    
    # Class multicfg
    class multicfg():

        def __init__(self):
            self.s2t = {}
            self.alltraits = []
            self.trait2fg = {}
            self.trait2bg = {}

        def update_dictionary(self, traitname, species, group):
            try:
                self.s2t[species].append(traitname + "_" + group)
            except:
                self.s2t[species] = [traitname + "_" + group]
                
            if group == "1":
                try:
                    self.trait2fg[traitname].append(species)
                except:
                    self.trait2fg[traitname] = [species]

            if group == "0":
                try:
                    self.trait2bg[traitname].append(species)
                except:
                    self.trait2bg[traitname] = [species]

    z = multicfg()

    # Step 1: read the trait config file and fetch fg and bg sizes
    with open(trait_config_file) as tcf_handle:
        tcf = tcf_handle.read().splitlines()
    
    fg = []
    bg = []

    for l in tcf:
        try:
            c = l.split("\t")
            if c[1] == "1":
                fg.append(c[0])
            elif c[1] == "0":
                bg.append(c[0])
        except:
            pass
    
    fg_len = len(fg)
    bg_len = len(bg)

    # Step 2: import the tree
    tree = dendropy.Tree.get(path=tree_file, schema=tree_format)

    # WORKFLOW 1: bootstrap in random mode

    # Extract the species

    if mode == "random":
        species = []
        for i, t1 in enumerate(tree.taxon_namespace[:-1]):
            species.append(t1.label.replace(" ", "_"))

        allcycles = []

        # Run the bootstrap
        for i in range(1,cycles):
            cycle_tag = "b_" + str(i)

            z.alltraits.append(cycle_tag)

            simulated_fg = []
            simulated_bg = []

            while len(simulated_fg) < fg_len:
                k = random.choice(species)

                if k not in simulated_fg:
                    simulated_fg.append(k)
        
            while len(simulated_bg) < bg_len:
                k = random.choice(species)
                if k not in simulated_bg and k not in simulated_fg:
                    simulated_bg.append(k)
            
            allcycles.append([simulated_fg, simulated_bg, cycle_tag])
            

        for c in allcycles:

            z.alltraits.append(c[2])
            # Foreground
            for s in c[0]:
                z.update_dictionary(c[2], s, "1")

            # Background
            for s in c[1]:
                z.update_dictionary(c[2], s, "0")
    return z
        
        


# Test
'''
tree_path = "bt_testing/traits.speciestree.210512.nh"
cfg_path = "tests/AOB_config.txt"

simtrait(cfg_path, tree_path)
'''