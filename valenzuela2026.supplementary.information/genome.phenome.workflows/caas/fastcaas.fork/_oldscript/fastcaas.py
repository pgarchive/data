#!/usr/bin/env python3

program_info = '''
FAST CAAS v. 0.7 (beta)                    
FAST Convergent Amino Acid Substitution detection tool
'''

# Authors:  Fabio Barteri   -    Universitat Pompeu Fabra (http://upf.edu)
#                                Institut de Biologia Evolutiva (https://www.ibe.upf-csic.es/)
#                                fabio.barteri@upf.edu
#           
#           Xavier Farré    -    Universitat Pompeu Fabra (http://upf.edu)
#                                Institut de Biologia Evolutiva (https://www.ibe.upf-csic.es/)

'''
TABLE OF CONTENTS
### OPTION PARSER AND MAIN SETTINGS
# Option parser
#   INPUTS
#       Run on single alignment
#       Run on an alignments folder
#       Config file
#   FILTERING
#       1. limit output to scenario (default: all scenarios)
#       2. filter per gaps 
#   OUTPUT
#   PARALLELISATION
# Checks whether the options are correct
# Settings (transferring parsed options to variables)
# Filtering settings (transferring parsed options to variables)
### FUNCTIONS
#   cfgparse() - Parses the config file
#   open_alignment() - Opens the alignment file
#   position_info() - information on the position
#   label_convergent() - labels the convergent substitutions
#   scan_alignment() - scans the single alignments
    # Step 1: filter the sequences in the cfg
    zeros_subset = list(filter(lambda x: x.id in cfgzeros, imported_alignment_file_ali))
    # Step 2:  scan the alisubset for mutations, gaps and all the rest
            # GAPS FILTERING
# chuncks() - divides the list in chuncks (n-sized shorter lists) for parallelisation
# chunk_scan() - scans the files in a single chunk
# print_header() - prints a header line on output file
### PRINT THE HEADER ON OUTPUT FILE (AND INIT IT)
### WORKFLOW ON A SINGLE ALIGNMENT
### WORKFLOW ON AN INPUT FOLDER (parallelised)                                                      
'''


### OPTION PARSER AND MAIN SETTINGS
############################################################

import sys,os
from optparse import OptionParser

# Option parser

parser = OptionParser()

#   INPUTS

#   Run on single alignment
parser.add_option("-a", "--alignment", dest="single_alignment",
                  help="Run caasper on a single alignment.", default = "none")
#   Run on an alignments folder
parser.add_option("-f", "--folder", dest="alignments_folder",
                  help="Run caasper on a folder containing more alignments. (-f/--folder $path_to_folder", default = "none")
#   Alignment format
parser.add_option("--fmt", dest="ali_format",
                  help="File format of the alignment (accepted formats in documentation). Default: fasta. Accepted: clustal, emboss, fasta, \
                  fasta-m10, ig, maf, mauve, msf, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm.", default = "fasta")
#   Config file
parser.add_option("-c", "--config", dest="config_file",
                  help="The config file (read documentation). (-c/--config $path_to_config_file", default = "none")

#   FILTERING

#       1. limit output to scenario (default: all scenarios)
parser.add_option("-s", "--scenarios", dest="scenarios_string",
                  help="Limit the result to some scenarios. Scenarios are indicated with numbers from 1 to 4 and must \
                      be provided as comma separated (e.g.: -s 1,2,3). See documentation for more details on the scenarios", default = "1,2,3,4")

#       2. filter per gaps 
parser.add_option("-g", "--maxgaps", dest="maxgaps_string",
                  help="Filter by gaps ratio (between 0 and 1). This option applies to both zero and uno groups \
                      e.g. -g 0.5 will accept positions with GapsRatio <= 0.5. Default = 1 (all positions accepted).", default = "1")

parser.add_option("-z", "--max_zero_gaps", dest="max_zero_gaps_string",
                  help="Filter by gaps ratio (between 0 and 1). This option applies to the zero group only\
                      e.g. -z 0.5 will accept positions with GapsRatio <= 0.5. Default = 1 (all positions accepted).", default = "1")

parser.add_option("-u", "--max_uno_gaps", dest="max_uno_gaps_string",
                  help="Filter by gaps ratio (between 0 and 1). This option applies to the uno group only\
                      e.g. -u 0.5 will accept positions with GapsRatio <= 0.5. Default = 1 (all positions accepted).", default = "1")



# OUTPUT

parser.add_option("-o", "--out", dest="output_file",
                  help="Path to output file.", default = "fastcaas.output.tab")

# PARALLELISATION

parser.add_option("-k", "--cpus", dest="number_of_cpus",
                  help="Number of CPUs. Applies to folder-based usage only.", default = "none")


parser.usage = '''\n[ON SINGLE ALIGNMENT FILE]\nfastcaas.py -a $input_alignment -c $config_file -o $output_file\n
[ON FOLDER WITH SEVERAL ALIGNMENT FILES]\nfastcaas.py -f $the_input_folder/ -c $config_file -o $output_file'''


(options, args) = parser.parse_args()

# Checks whether the options are correct

if len(sys.argv) == 1:
    print(program_info)
    parser.print_help()
    exit()

if options.single_alignment == "none" and options.alignments_folder == "none":
    print(program_info)
    print("\n\n****ERROR: no input alignment(s) provided****\n\n")
    parser.print_help()
    exit()

elif options.single_alignment != "none" and options.alignments_folder != "none":
    print(program_info)
    print("\n\n****ERROR: you must choose either a single alignment (-a) or to loop on a whole folder (-f) ****\n\n")
    parser.print_help()
    exit()

if options.config_file == "none":
    print(program_info)
    print("\n\n****ERROR: no config file provided****\n\n")
    parser.print_help()
    exit()


# Settings (transferring parsed options to variables)

alignment = options.single_alignment
alignments_folder = options.alignments_folder
alignments_format = options.ali_format
cfgfile = options.config_file
output_file = options.output_file

cfgtag = cfgfile.split("/")[-1]

# Filtering settings (transferring parsed options to variables)

admitted_scenarios = []
ascen = options.scenarios_string.split(",")
for x in ascen:
    admitted_scenarios.append("scenario" + x)


maxgaps = float(options.maxgaps_string)
max_zero_gaps = float(options.max_zero_gaps_string)
max_uno_gaps = float(options.max_uno_gaps_string)


### FUNCTIONS
############################################################

from Bio import AlignIO
import multiprocessing as mp

#   cfgparse() - Parses the config file

def cfgparse(infile):
    class gspec():
        def __init__(self):
            self.zeros = []
            self.unos = []
    
    z = gspec()
    with open(infile, "r") as handle_file:
        handle = handle_file.read().split("\n")
    
    for line in handle:
        try:
            c = line.split("\t")
            if c[1] == "0":
                z.zeros.append(c[0])
            elif c[1] == "1":
                z.unos.append(c[0])
        except:
            pass
    return z

#   open_alignment() - Opens the alignment file
def open_alignment(alignment_file):

    class processed_alignment():
        def __init__(self):
            self.ali= []
            self.len = 0
            self.name = ""

    z = processed_alignment()

    try:
        z.ali = list(AlignIO.read(alignment_file, alignments_format))
        z.len = AlignIO.read(alignment_file, alignments_format).get_alignment_length()
    except:
        print("ERROR:", alignment_file, "was not recognised. You need to specify the alignment format by --fmt (clustal, fasta, phylip...). Please read documentation for supported alignment formats.")

    if "/" in alignment_file:
        z.name = alignment_file.split("/")[-1]
    else:
        z.name = alignment_file
    return z


#   position_info() - information on the position
def position_info(position_list):

    class posinfo():
        def __init__(self):
            self.tag = ""
            self.list = []
            self.ungapped = []
            self.ungapped_tag = ''
            self.gaps = 0

    zz = posinfo()
    for x in position_list:
        if x[1] != "-":
            zz.ungapped.append(x[0])
            if x[1] not in zz.list:
                zz.list.append(x[1])
        else:
            zz.gaps = zz.gaps +1

    
    zz.list.sort()
    zz.tag = "/".join(zz.list)
    zz.ungapped.sort()
    zz.ungapped_tag = ",".join(zz.ungapped)
    return zz

#   label_convergent() - labels the convergent substitutions
def label_convergent(list1, list2, thescenarios):

    class thelabel():
        def __init__(self):
            self.flag = True
            self.scenario = ""

    l = thelabel()

    for x in list1:
        if x in list2:
            l.flag = False
            break

    if len(list1) < 1 or len(list2) < 1:
        l.flag = False
    
    if len(list1) == 1 and len(list2) == 1:
        l.scenario = "scenario1"
    elif len(list1) == 1 and len(list2) > 1:
        l.scenario = "scenario2"
    elif len(list1) > 1 and len(list2) == 1:
        l.scenario = "scenario3"
    elif len(list1) > 1 and len(list2) > 1:
        l.scenario = "scenario4"

    if l.scenario not in thescenarios:
        l.flag = False
    
    return l


#   scan_alignment() - scans the single alignments
def scan_alignment(imported_alignment_file_ali, imported_alignment_file_len, imported_alignment_file_name, cfgzeros, cfgunos, thecfgtag):

    class detected_position():
        def __init__(self):
            self.alignment_name = ""
 
            self.position = ""
            self.zeros = ""
            self.unos = ""
            self.scenario = ""

            self.zero_species_in_config = ""
            self.uno_species_in_config = ""

            self.zero_found_in_alignment = ""
            self.uno_found_in_alignment = ""

            self.zero_ungapped_species_number = ""
            self.uno_ungapped_species_number = ""

            self.zeros_gaps = ""
            self.unos_gaps = ""

            self.zeros_gaps_ratio = ""
            self.unos_gaps_ratio = ""

            self.total_gaps = ""
            self.total_gaps_ratio = ""

            self.zeros_species = ""
            self.unos_species = ""


        
        def return_line(self):

            outline = "\t".join([
                thecfgtag,
                self.alignment_name,

                self.position,
                self.zeros,
                self.unos,
                self.scenario,

                self.zero_species_in_config,
                self.uno_species_in_config,

                self.zero_found_in_alignment,
                self.uno_found_in_alignment,

                self.zero_ungapped_species_number,
                self.uno_ungapped_species_number,

                self.zeros_gaps,
                self.unos_gaps,

                self.zeros_gaps_ratio,
                self.unos_gaps_ratio,

                self.total_gaps,
                self.total_gaps_ratio,

                self.zeros_species,
                self.unos_species,



            ])

            return outline

    # Step 1: filter the sequences in the cfg
    zeros_subset = list(filter(lambda x: x.id in cfgzeros, imported_alignment_file_ali))            ### Filter() function instead of for cycle
    unos_subset = list(filter(lambda x: x.id in cfgunos, imported_alignment_file_ali))              
    
    # Step 2:  scan the alisubset for mutations, gaps and all the rest

    output_lines = []

    for i in range(0, imported_alignment_file_len):
        zeros_subset_info = list(map(lambda x : (x.id, x.seq[i]), zeros_subset))
        unos_subset_info = list(map(lambda x : (x.id, x.seq[i]), unos_subset))

        zeros_position_info = position_info(zeros_subset_info)
        unos_position_info = position_info(unos_subset_info)

        flag = label_convergent(zeros_position_info.list, unos_position_info.list, admitted_scenarios)

        if flag.flag == True:
            strike = detected_position()

            strike.alignment_name = imported_alignment_file_name
            strike.position = str(i + 1)
            strike.zeros = zeros_position_info.tag
            strike.unos = unos_position_info.tag
            strike.scenario = flag.scenario

            strike.zero_species_in_config = str(len(cfgzeros))
            strike.uno_species_in_config = str(len(cfgunos))

            strike.zero_found_in_alignment = str(len(zeros_subset))
            strike.uno_found_in_alignment = str(len(unos_subset))

            strike.zero_ungapped_species_number = str(len(zeros_position_info.ungapped_tag.split(",")))
            strike.uno_ungapped_species_number = str(len(unos_position_info.ungapped_tag.split(",")))


            strike.zeros_species = zeros_position_info.ungapped_tag
            strike.unos_species = unos_position_info.ungapped_tag

            strike.zeros_gaps = str(zeros_position_info.gaps)
            strike.unos_gaps = str(unos_position_info.gaps)

            strike.zeros_gaps_ratio = str(zeros_position_info.gaps / float(len(zeros_subset)))
            strike.unos_gaps_ratio = str(unos_position_info.gaps / float(len(unos_subset)))

            totalgaps = zeros_position_info.gaps + unos_position_info.gaps
            totalspecies = len(cfg.zeros) + len(cfg.unos)

            strike.total_gaps = str(totalgaps)
            strike.total_gaps_ratio = str(totalgaps / float(totalspecies))

            # GAPS FILTERING

            gapsflag = True

            if totalgaps / float(totalspecies) > maxgaps:
                gapsflag = False

            if zeros_position_info.gaps / float(len(cfg.zeros)) > max_zero_gaps:
                gapsflag = False

            if unos_position_info.gaps / float(len(cfg.unos)) > max_uno_gaps:
                gapsflag = False

            if len(strike.return_line()) > 1 and gapsflag == True:
                output_lines.append(strike.return_line())
    
    return "\n".join(output_lines)


# chuncks() - divides the list in chuncks (n-sized shorter lists) for parallelisation

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

# chunk_scan() - scans the files in a single chunk

def chunk_scan(single_chunk, cfgzeros, cfgunos):
    outlist = []
    for file in single_chunk:
        f = open_alignment(file)
        outlist = outlist + scan_alignment(f.ali, f.len, f.name, cfgzeros, cfgunos, cfgtag)
    return("\n".join(outlist))

# print_header() - prints a header line on output file
# [STILL TO WRITE]

def print_header(output_file):

    theout = open(output_file, "w")

    outline = "\t".join([
        "ConfigName",
        "AlignmentName",
        "Position",
        "0_GroupAA",
        "1_GroupAA",
        "Scenario",
        "0_species_in_config",
        "1_species_in_config",
        "0_found_in_alignment",
        "1_found_in_alignment",
        "0_ungapped_species_n",
        "1_ungapped_species_n",
        "0_NumberOfGaps",
        "1_NumberOfGaps",
        "0_GapsRatio",
        "1_GapsRatio",
        "TotalGaps",
        "Total_GapsRatio",
        "0_ungapped_species",
        "1_ungapped_species",

    ])

    print(outline, file=theout)
    theout.close()

### PRINT THE HEADER ON OUTPUT FILE (AND INIT IT)
############################################################

print_header(output_file)

### WORKFLOW ON A SINGLE ALIGNMENT
############################################################
if alignment != "none":
    print(program_info)
    print("Processing", alignment)
    out = open(output_file, "a")
    cfg = cfgparse(cfgfile)
    ali = open_alignment(alignment)
    outlines = scan_alignment(ali.ali, ali.len, ali.name, cfg.zeros, cfg.unos, cfgtag)
    out.write(outlines)
    out.close()
    print("DONE!\nResults are available at", output_file, "\n\n")

### WORKFLOW ON AN INPUT FOLDER (parallelised)                                                      ### PARALLELISED WORKFLOW ON n CPUs
############################################################

elif alignments_folder != "none":

    import glob
    import multiprocessing as mp


    # Setting up
    print(program_info)
    print("Processing", alignments_folder)
    cfg = cfgparse(cfgfile)
    thealignments = glob.glob(alignments_folder + "/*".replace("//","/"))
    imported_alginments = list(map(open_alignment, thealignments))

    # Unisci open alignments e process alignments e poi parallelizza a partire dai chunks

    if options.number_of_cpus == "none":
        pool = mp.Pool(mp.cpu_count())
    else:
        try:
            pool = mp.Pool(int(options.number_of_cpus))
            if int(options.number_of_cpus) > mp.cpu_count():
                print("WARNING. TOO MANY CPUs PROVIDED. FAST CAAS WILL RUN ON ALL AVAILABLE CPUs")
        except:
            print("WARNING. NON READABLE NUMBER OF CPUS PROVIDED. FAST CAAS WILL RUN ON ALL AVAILABLE CPUs")
            pool = mp.Pool(mp.cpu_count())

    results = [pool.apply(scan_alignment, args=(alignment.ali, alignment.len, alignment.name,cfg.zeros, cfg.unos,cfgtag)) for alignment in imported_alginments]
    pool.close()

    out = open(output_file, "a")
    for x in results:
        if len(x) > 2:
            print(x, file= out)
    out.close()

print("Done. Results are available in", output_file)