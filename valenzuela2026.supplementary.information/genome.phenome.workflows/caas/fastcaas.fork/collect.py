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

MODULE NAME: collect.py
DESCRIPTION: contains the two functions to collect and filter the results
DEPENDENCIES: none 
CALLED BY: __main__ 
'''

import glob,functools

###  Function filter_results. Filters the results from several CAAS files, given a max pvalue (default = 0.05).

def filter_single_result(resultfile, maxpval):

    def isok(string, max_pval):
        out = False
        try:
            c = string.split("\t")
            pv = float(c[4])

            if pv <= maxpval:
                out = True
        except:
            pass

        return out

    with open(resultfile) as rf_handle:
        rf = rf_handle.read().splitlines()
    
    valids = filter(functools.partial(isok, max_pval = maxpval), rf)
    return "\n".join(list(valids))

### FUNCTION filter_multiple_results(). Filters more than one result

def filter_multiple_results(results_folder, the_max_pval, output_file):
    filelist = glob.glob(results_folder + "/*caas")

    valid_results = map(functools.partial(filter_single_result, maxpval = the_max_pval), filelist)
    valid_results = filter( lambda x : len(x) > 2, valid_results)
    k = list(valid_results)
    
    of = open(output_file, "w")
    print("\n".join(k), file = of)
    of.close()



def mergeboot(results_folder, output_folder, pvalue):

    class excess_per_cfg():

        def __init__(self):
            self.cfgnames = []
            self.name2all = {}
            self.name2sig = {}
        
        def import_boot_results(self, results_folder, pvalue):

            ali_boots = glob.glob(results_folder + "/results_per_alignment/*boot")

            for a in ali_boots:

                with open(a) as handle:
                    imported_ali = handle.read().splitlines()
                
                for line in handle:
                    try:
                        c = line.split("\t")

                        pv = float(c[2])
                        pos = c[0]
                        cfgname = c[4].split("/")[-1]

                        # Add cfgname to the list
                        if cfgname not in self.cfgnames:
                            self.cfgnames.append(cfgname)
                        
                        # Compare the pvalues and add the position

                        if pv <= pvalue:

                            try:
                                self.name2sig[cfgname].append(pos)
                            except:
                                self.name2sig[cfgname] = [pos]

                            try:
                                self.name2all[cfgname].append(pos)
                            except:
                                self.name2all[cfgname] = [pos]
                        else:

                            try:
                                self.name2all[cfgname].append(pos)
                            except:
                                self.name2all[cfgname] = [pos]

                    except:
                        pass

        
        def printresume(self, output_folder):

            # A clean option
            for x in self.cfgnames:
                out = open()


    
    # Select the alignments
    ali_boots = glob.glob(results_folder + "/results_per_alignment/*boot")
    

    
    return 0

### FUNCTION croupier(). Updates a dictionary to sort all the mutations per each bootstrap

def croupier(infile, pvalue, dictionary):
    # Printing info
    print("merging boostrap for", infile.split("/")[-1].split(".")[0])
    
    with open(infile) as h:
        f = h.read().splitlines()
    outtag = ""
    for l in f:
        try:
            c = l.split("\t")
            pv = float(c[2])
            pos = c[0]
            outtag = c[4].split("/")[-1]

            if pv <= pvalue:

                positives = c[3].split(",")
                for p in positives:
                    if len(p) > 1:
                        try:
                            dictionary[p].append(pos)
                        except:
                            dictionary[p] = [pos]
                            pass
        except:
            pass
    return outtag

### Function merge_bootstrap()

def merge_bootstrap(working_folder, pvalue):

    allboots = {}
    underpval = {}


    # Read the bootstrap values

    ali_boots = glob.glob(working_folder + "/results_per_alignment/*boot")

    for b in ali_boots:
        w = croupier(b, pvalue, underpval)
        x = croupier(b, 1, allboots)
    
    # Print the bootstrap resume

    outfile = open(working_folder + "/" + w + "__bootstrap.table", "w")
    
    print("\t".join(["Cycle", "All", "UnderPvalue_" + str(pvalue), "All_mutations", "UnderPvalue_" + str(pvalue) + "_mutations"]), file=outfile)

    for k in allboots.keys():
        try:
            upv = underpval[k]
        except:
            upv = []
        
        print(
            "\t".join([k, str(len(allboots[k])), str(len(upv)), ",".join(allboots[k]), ",".join(upv)]),
            file = outfile
        )
    outfile.close()

### TEST
'''
f = "AOB_config.txt_caas/results_per_alignment/ZACN.Homo_sapiens.filter2.fasta.caas"
d = "AOB_config.txt_caas/results_per_alignment/"
pval = 0.05

#merge_bootstrap("AOB_config.txt_caas/", pval)
filter_multiple_results(d, pval, output_file="ghghgh.caas.tuf")
'''