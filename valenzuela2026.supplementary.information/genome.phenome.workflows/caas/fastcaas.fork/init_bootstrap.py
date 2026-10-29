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

from warnings import simplefilter
import random, os, sys

# FUNCTION simtrait() Simulate trait function

def simtrait(trait_config_file, species_file, mode = "random", cycles = 1000, traits_outfile = "none", tree_file = "none", phenotype_values_file = "none", permulation_selection_strategy = "random"):
    
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
        
        def print_traits(self, outfile):

            o = open(outfile, "w")
            for x in self.trait2bg.keys():
                print("\t".join([   x,
                                    ",".join(self.trait2fg[x]),
                                    ",".join(self.trait2bg[x])
                                    ]), file = o)
            o.close()

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

    # Step 2: import the species

    with open(species_file) as species_handle:
        species_lines = species_handle.read().splitlines()
    
    # Step 3: import the species with a phenotype

    with open(phenotype_values_file) as phenotypes_handle:
        phenotypes_lines = phenotypes_handle.read().splitlines()

        species_pool = []
        for line in phenotypes_lines:
            try:
                c = line.split("\t")
                species_pool.append(c[0])
            except:
                pass

    # WORKFLOW 1: bootstrap in random mode


    if mode == "random":
        
        # Extract the species

        species = []

        for line in species_lines:
            try:
                c = line.split("\t")
                if c[0] in species_pool:
                    species.append(c[0])

            except:
                pass

        allcycles = []

        # Run the bootstrap
        for i in range(1,cycles+1):
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
    
    elif mode == "phylogeny-guided":

        # Extract the species dictionaries

        s2g = {}
        g2s = {}

        for line in species_lines:
            try:
                c = line.split("\t")
                if c[0] in species_pool:
                    # Update s2g
                    s2g[c[0]] = c[1]

                    # Update g2s
                    try:
                        g2s[c[1]].append(c[0])
                    except:
                        g2s[c[1]] = [c[0]]
            except:
                pass
        
        species = s2g.keys()

        allcycles = []

        # Run the bootstrap
        for i in range(1,cycles+1):
            cycle_tag = "b_" + str(i)

            z.alltraits.append(cycle_tag)

            simulated_fg = []
            simulated_bg = []

            for x in fg:
                flag = 0
                while flag != 1:
                    g = s2g[x]
                    x = random.choice(g2s[g])

                    if x not in simulated_fg:
                        simulated_fg.append(x)
                        flag = 1


            for x in bg:
                flag = 0
                while flag != 1:
                    g = s2g[x]
                    x = random.choice(g2s[g])

                    if x not in simulated_bg and x not in simulated_fg:
                        simulated_bg.append(x)
                        flag = 1


            simulated_bg.sort()
            simulated_fg.sort()
            
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
    
    elif mode == "permulations":

        if tree_file == "none":
            print("ERROR: for permulations mode you must provide a rooted and fully dichotomic tree file in newick format.")
            print("Exiting.")
            exit()
        if phenotype_values_file == "none":
            print("ERROR: for permulations mode you must provide a tsv file with phenotype values (1/0 for binary traits).")
            print("Exiting.")
            exit()
        
        script_path = os.path.realpath(sys.argv[0]).split("/")
        script_path.remove(script_path[-1])
        rscript_path = "/".join(script_path) + "/permulations.r"

        os.system("ln -s " + rscript_path)
        
        r_line = " ".join([
            
            "./permulations.r",                     # Rscript
            tree_file,                              # args[1]
            trait_config_file,                      # args[2]
            str(cycles),                            # args[3]
            permulation_selection_strategy,         # args[4]
            phenotype_values_file,                  # args[5]
            traits_outfile                          # args[6]
        ])
        os.system(r_line)

# FUNCTION simtrait_revive() revive simulated trait from

def simtrait_revive(traitfile):
    
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
            
            self.alltraits.append(traitname)
        
        
        def print_traits(self, outfile):
            o = open(outfile, "w")
            for x in self.trait2bg.keys():
                print("\t".join([   x,
                                    ",".join(self.trait2fg[x]),
                                    ",".join(self.trait2bg[x])
                                    ]), file = o)
            o.close()

        
    # Declare multicfg instance

    z = multicfg()

    # Open the traitfile

    with open(traitfile) as tf_handle:
        tf = tf_handle.read().splitlines()
    
    for line in tf:
        try:
            c = line.split("\t")
            cycleid = c[0]
            fg = c[1].split(",")
            bg = c[2].split(",")

            # Foreground update

            for s in fg:
                z.update_dictionary(cycleid, s, "1")


            # Background update

            for s in bg:
                z.update_dictionary(cycleid, s, "0")

        except:
            pass

    return z



'''
# Test
species_path = "_tests/sp2fam.210727.tab"
cfg_path = "_tests/AOB_config.txt"

st = simtrait(cfg_path, species_path, mode="phylogeny-guided")
st.print_traits("provaoutfile.tuf")

rst = simtrait_revive("provaoutfile.tuf")
rst.print_traits("provaoutfile.revived.tuf")
'''