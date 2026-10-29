#!/usr/bin/env python3

# Fastcaas - fetchmismatches.py module.
# This module creates an index out of a group of alignments.

'''
Info.
Fetchmismatches is deployed to identify pairwise mismatches between species out of a MSA. The MSA
is scanned.
'''

import glob, functools

### SETTINGS

configs_folder = "config_files.210528"

### FUNCTIONS

# FUNCTION update dictionary

def update_dictionary(dictionary, key, value):
    try:
        dictionary[key].append(value)
    except:
        dictionary[key] = [value]

# FUNCTION load multi cfg dictionary
# Loads the multi cfg dictionary


def load_cfg_dictionary(dictionary_folder):

    class multicfg():

        def __init__(self):
            self.s2t = {}
            self.alltraits = []

        def update_dictionary(self, key, value):
            try:
                self.s2t[key].append(value)
            except:
                self.s2t[key] = [value]

    z = multicfg()

    for x in glob.glob(dictionary_folder + "/*"):

        traitname = x.split("/")[-1]
        z.alltraits.append(traitname)

        with open(x) as singlecfg_f:
            singlecfg =  singlecfg_f.read().splitlines()
        
        for line in singlecfg:
            try:
                c = line.split()
                z.update_dictionary(c[0], traitname + "_" + c[1])
            except:
                pass

    return z