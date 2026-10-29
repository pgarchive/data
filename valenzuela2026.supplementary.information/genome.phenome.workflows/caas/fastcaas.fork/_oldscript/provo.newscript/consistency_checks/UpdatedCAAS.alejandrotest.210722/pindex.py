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

MODULE NAME: pindex.py
DESCRIPTION: phenotype indexing of a single trait or multiple traits
INPUTS:      single or multiple trait binary file

TABLE OF CONTENTS
------------------------------------------
update_dictionary()         a function to update a dictionary with
                            new information. No, there is no built-in method for
                            this.

load_cfg_dictionary()       Loads the multi cfg dictionary

'''

import glob

# FUNCTION update dictionary
# A function to update a dictionary with new information. No, there is no built-in method for this.

def update_dictionary(dictionary, key, value):
    try:
        dictionary[key].append(value)
    except:
        dictionary[key] = [value]

# FUNCTION load multi cfg dictionary
# Loads the multi cfg dictionary

def load_cfg(input_path, mode = "mono"):

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

    if mode == "multi":

        for x in glob.glob(input_path + "/*"):

            traitname = x.split("/")[-1]
            z.alltraits.append(traitname)

            with open(x) as singlecfg_f:
                singlecfg =  singlecfg_f.read().splitlines()
            
            for line in singlecfg:
                try:
                    c = line.split()
                    z.update_dictionary(traitname, c[0], c[1])
                except:
                    pass
    
    elif mode == "mono":

        traitname = input_path.split("/")[-1]
        z.alltraits.append(traitname)

        with open(input_path) as singlecfg_f:
            singlecfg =  singlecfg_f.read().splitlines()
            
        for line in singlecfg:
            try:
                c = line.split()
                z.update_dictionary(traitname, c[0], c[1])
            except:
                pass

    return z