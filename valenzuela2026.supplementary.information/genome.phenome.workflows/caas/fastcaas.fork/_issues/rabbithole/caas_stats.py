#!/usr/bin/env python3

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

MODULE NAME: caas_stats.py
DESCRIPTION: Hypergeometric pvalue assignation.
INPUTS:      Input MSAs
CALLED BY:   caas_id.py, fastcaas_core.py


THIS IS THE RABBIT HOLE 


'''                                                       
import time

st_imports = time.time()

from itertools import combinations
import random
import functools
from scipy import stats as ss
import glob

et_imports = time.time() - st_imports

print("Done. " + str(et_imports) + " seconds elapsed")

# FUNCTION count_symbols() - Counts the symbols (AAs) 
def count_symbols(list1, list2):
    outdict = {}
    for x in list1:
        outdict[x] = list2.count(x)
    
    return outdict

# FUNCTION pstate() - Calculate the statistics of a symbol combination ()

def pstate(iterable, freq_dictionary, set_size, contrast_size):

    '''
    Hypergeometric variables

    M = IL MAZZO DI CARTE
    Mn = IL MAZZO DI CARTE SAPENDO CHE L'AVVERSARIO (BACKGROUND) HA COLTO LE SUE (M - CONTRAST_SIZE)
    N = IL NUMERO DI CARTE CHE VUOI NEL MAZZO = sum of the frequences of the iterable
    n = LE CARTE CHE HAI IN MANO (SET_SIZE)
    k LE CARTE CHE HAI IN MANO E CHE COINCIDONO CON (N) -----> n = k (SET_SIZE)


    '''
    
    M = functools.reduce(lambda a, b: a+b, freq_dictionary.values())
    Mn = M - contrast_size

    it_freqs = map(lambda x : freq_dictionary[x], iterable)
    N = functools.reduce(lambda a, b: a+b, it_freqs)
    n = set_size
    k = set_size

    if N > Mn:
        p = 0
    else:
        l = ss.hypergeom(Mn, n, N)
        p = l.pmf(k)

    return p

### FUNCTION sstate() - Returns

def sstate(combination, freq_dictionary, fg_size, bg_size):
    f = pstate(list(combination[0]), freq_dictionary, fg_size, bg_size)
    b = pstate(list(combination[1]), freq_dictionary, bg_size, fg_size)

    return f * b

# FUNCTION calcpval_random() - fetches a pvalue from a line, given foreground and background sizes
def calcpval_random(line_dictionary, genename, fg_size, bg_size):

    class lstats():
        def __init__(self):
            self.p = ""

            # Symbols
            self.s =[]                      # Symbols. All the amino acids in the position.
            self.sfreq = {}                 # Frequency per each symbol in th position
            self.c = []                     # One to many combinations
            self.uc = []                    # One to one combinations

            self.fgsize = fg_size
            self.bgsize = bg_size

            self.pvalue = 0

            self.out_tuple = ()

        def combine(self):
            for x in self.s:
                bg = list(set(self.s))
                bg.remove(x)
                self.c.append([[x], bg])
                self.c.append([bg, [x]])
            
            self.uc = list(combinations(self.s, 2))
            rec = map(lambda x : (x[1], x[0]), self.uc)

            self.uc = self.uc + list(rec)
        
        def dostat(self):

            combresults = list(map(functools.partial(sstate, freq_dictionary = self.sfreq, fg_size = self.fgsize, bg_size = self.bgsize), self.c))
            p_comb = functools.reduce(lambda a, b: a+b, combresults)
                        
            if len(self.s) > 2:                     # You will subtract one to one combinations only w 3+ AAs
                
                u_combresults = list(map(functools.partial(sstate, freq_dictionary = self.sfreq, fg_size = self.fgsize, bg_size = self.bgsize), self.uc))
                up_comb = functools.reduce(lambda a, b: a+b, u_combresults)
                self.pvalue = p_comb - up_comb

            else:
                self.pvalue = p_comb
            
            self.out_tuple = (genename + "@" + str(self.p), self.pvalue)

    z = lstats()

    symbols = list(map(lambda x : x.split("@")[0], line_dictionary.values()))
    z.p = list(map(lambda x : x.split("@")[1], line_dictionary.values()))[0]
    
    try:
        symbols.remove("-") # Clean the gaps...
    except:
        pass                # ...if any.

    # Load the z object

    z.s = set(symbols)
    z.sfreq = count_symbols(z.s, symbols)
    z.combine()
    z.dostat()

    return z.out_tuple

### FUNCTION genepval() - iterates the pvalue calculation over a whole gene

def genepval(imported_alignment, gpv_fg_size, gpv_bg_size, mode = "random"):
    pval_dictionary = {}

    if mode == "random":
        pvals = map(functools.partial(calcpval_random, genename = imported_alignment.genename, fg_size = gpv_fg_size, bg_size = gpv_bg_size), imported_alignment.d)
        pval_dictionary = dict(pvals)

    return pval_dictionary


### FUNCTION sim_bootstrap() - simulates a bootstrap starting from the generated pvalues

def sim_bootstrap(pvaldictionary, number_of_cycles = 1000):
    # Step 1: Create the bootstrap values
    bvals = list(map(lambda x : "b_" + str(x), (x for x in range(1,number_of_cycles))))

    # Step 2: Calculate the number of successes per mutation
    nfold_map = map(lambda x : (x[0], int(x[1] * number_of_cycles)), pvaldictionary.items())
    successes = dict(nfold_map)

    # Step 3: scoop a random value
    xtraction = map(lambda x: random.sample(bvals, x), successes.values())
    xtr_results = functools.reduce(lambda a, b: a+b, xtraction)
    single_count = map(lambda x: (x, xtr_results.count(x)), bvals)
    bootstrap_dictionary = dict(single_count)

    return bootstrap_dictionary




### Testing 3: The craziest idea

from alimport import *

# Step 1. Import the simulated values
with open("test.littersize.107genes/simulated.results.trigly.tab") as sim_h:
    sims = sim_h.read().splitlines()

sims_t = map(lambda x : x.split("\t"), sims)
pos2spv = dict(sims_t)

# Step 2. calculate the values for the genes (limited to 5) for dev reasons.

# The FG and BG sizes based on Triglycerides

foreground_size = 2
background_size = 6

slice_t = min(background_size, foreground_size)

random_alignments = glob.glob("test.littersize.107genes/genes/*fasta")

pos2cpv = {}

for x in random_alignments:
    st_stat= time.time()
    
    slicex = slice(x, column_threshold=slice_t)
    slicex_pvals = genepval(slicex, foreground_size, background_size)
    pos2cpv.update(slicex_pvals)    
    et_stat = time.time() - st_stat
    print(x + " Done. " + str(len(slicex.d))  + " Analysed in " + str(et_stat) + " " + str(len(slicex.d)/et_stat) + " lines per second")

# Print the output

o = open("simcalc.comparison.trigly.tab", "w")

print("\t".join(["position","pvalue", "strategy"]), file = o)

for p in pos2cpv.keys():
    try:
        outline_1 = "\t".join([p, str(pos2cpv[p]), "Calculated"])
        outline_2 = "\t".join([p, str(pos2spv[p]), "Simulated"])
        print(outline_1, file=o)
        print(outline_2, file=o)
    except:
        pass

o.close()



'''
### TESTING 2:  the bootstrap simulation
# Data 1: the simulated bootstrap values

with open("trigly.boot.dict") as sim_bs_h:
    simdict = sim_bs_h.read().splitlines()

print("Importing the bootstrap from simulated values...")
sim_bs_map = map(lambda x : x.split("\t"), simdict)
boot2sim = dict (sim_bs_map)
# Data 2: the calcluated pvalue dictionary (NO WEIGHTS)

#   revive the pvalue dictionary

print("Importing the bootstrap from calculated values - NO WEIGHTS...")

with open("simcalc.comparison.trigly.tab") as sc_handle:
    sc = sc_handle.read().splitlines()

filtered_sc = filter(lambda x: x.split("\t")[-1] == "Calculated", sc)
mapped_sc = map(lambda x: (x.split("\t")[0], float(x.split("\t")[1])), filtered_sc)
pval_dict = dict(mapped_sc)


boot2calc = sim_bootstrap(pval_dict)

### TESTING 1: Building the calcpval_random() function

This testing has been made on 211027 on LitterSz_config.txt on 107 genes from genes.list file (in rabbithole folder).

Simulated values: pvalues from bootstrap on simulated traits
Calculated values: pvalues from hypergeometric calculation

from alimport import *

# Step 1. Import the simulated values
with open("test.littersize.107genes/simulated.results.trigly.tab") as sim_h:
    sims = sim_h.read().splitlines()

sims_t = map(lambda x : x.split("\t"), sims)
pos2spv = dict(sims_t)

# Step 2. calculate the values for the genes (limited to 5) for dev reasons.

# The FG and BG sizes based on Triglycerides

foreground_size = 2
background_size = 6

slice_t = min(background_size, foreground_size)

random_alignments = glob.glob("test.littersize.107genes/genes/*fasta")

pos2cpv = {}

for x in random_alignments:
    st_stat= time.time()
    
    slicex = slice(x, column_threshold=slice_t)
    slicex_pvals = genepval(slicex, foreground_size, background_size)
    pos2cpv.update(slicex_pvals)    
    et_stat = time.time() - st_stat
    print(x + " Done. " + str(len(slicex.d))  + " Analysed in " + str(et_stat) + " " + str(len(slicex.d)/et_stat) + " lines per second")

# Print the output

o = open("simcalc.comparison.trigly.tab", "w")

print("\t".join(["position","pvalue", "strategy"]), file = o)

for p in pos2cpv.keys():
    try:
        outline_1 = "\t".join([p, str(pos2cpv[p]), "Calculated"])
        outline_2 = "\t".join([p, str(pos2spv[p]), "Simulated"])
        print(outline_1, file=o)
        print(outline_2, file=o)
    except:
        pass

o.close()
'''