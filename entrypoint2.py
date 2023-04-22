
import json
import sys
import os
from locale import atof, setlocale, LC_NUMERIC
from datetime import datetime
import hashlib
import logging


from pathlib import Path

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import miRNA
# '''__all__ = []
# __version__ = 0.1
# __date__ = '2023-04-11'
# __updated__ = '2022-04-11'
#
# import sequence
# import miRNA
#
# DEBUG = 1
# TESTRUN = 0
# PROFILE = 0'''

def readFastaFile(filename):
    '''
    load specified fasta file and store header and sequence as entries in two lists
    :param self:
    :return:
    '''

    print("load sequences from fasta file <" + filename + ">")
    global headerLines
    global sequenceLines
    global speciesCode
  
  # List of species codes to process
    speciesCode = ("hsa", "mmu", "rno")

# Loop through each species code
    for code in speciesCode:
        # Set the input and output file paths
        input_file = f"input_{code}.fa"
        output_file = f"output_{code}.fa"
        
        # Run the command to process the input file
        #cmd = f"python process.py {C:\Users\user\AdvancedPrograming\day4\data\mature.fa} {output_file}"
    input_file = r"C:\Users\user\AdvancedPrograming\day4\data\mature.fa"
    output_file = f"output_{code}.fa"
    #cmd = f"python process.py {input_file} {output_file}"
    #os.system(cmd)

    # load the fasta lines into a list
    try:
        fFA = open(filename, 'r')
        fastaLines = fFA.readlines()
        fFA.close()
    except Exception as e:
        raise(e)

    headerLines = []
    headerLine = ""
    sequenceLines = []
    sequence = ""

    s = 0
    seqCount = 0
    for fastaLine in fastaLines:
        if fastaLine[0] == '>':
            seqCount +=1
            if s > 0 and headerLine.startswith(speciesCode):
                headerLines.append(headerLine)
                sequenceLines.append(sequence)
                sequence = ""
            headerLine = fastaLine[1:].strip()
            sequence = ""
            
        else:
            sequence = sequence + fastaLine.strip()
        s += 1
    if headerLine.startswith(speciesCode):
        headerLines.append(headerLine)
        sequenceLines.append(sequence)   
 
    # print("loaded <" + str(seqCount) + "> sequences and kept <"+ str(len(headerLines)) + "> with species code [" + speciesCode + "]")
    print("loaded <" + str(seqCount) + "> sequences and kept <"+ str(len(headerLines)) + "> with species code [" + str(speciesCode) + "]")
        
    return len(headerLines)

def getUniqueSeedSequences():
    '''
    get the unique seed sequences in the list of sequences loaded from the fasta file
    :return:
    '''
    global seedBegin
    global seedEnd
    
    print("get unique seed sequences from sequence list")
    print("seed region is defined to run from <" + str(seedBegin) + ">--><" + str(seedEnd) + ">")
    
    uniqSeedSeqs = []
    seqNo = 0
    for seqLine in sequenceLines:
        miRSeq = miRNA.MiRNA(headerLines[seqNo], seqLine)
        thisSeedSeq = miRSeq.getSeedSequence(seedBegin, seedEnd)

        if thisSeedSeq not in uniqSeedSeqs:
            uniqSeedSeqs.append(thisSeedSeq)


    print("found <" + str(len(uniqSeedSeqs))+ "> unique seed sequences")
    return uniqSeedSeqs


def writeUniqSeqs(uSeqs):
    '''
    write unique seed regions to an output file
    '''
    print("write unique seed sequences to fasta file")
    import os
    from pathlib import Path
    
    foldername = os.path.dirname(filename)
    basename = Path(filename).stem    
    outputfafile = os.path.join(foldername, basename + "__uniqseeds" + ".fa")
    
    print("output fasta file is <" + outputfafile + ">")
    file = open(outputfafile,'w')
    s = 1
    for uSeq in uSeqs:
        file.writelines(">uniqseed_" + str(s) + "\n")
        file.writelines(uSeq + "\n")
        s+=1
    file.close()    
    
    return uSeqs  

def getNucleotideFrequencyMatrix(uniqSeedSeqs):
    '''
    calculate nucleotide frequencies at each position
    '''
    if uniqSeedSeqs is None:
        raise ValueError("Input 'uniqSeedSeqs' is None.")
    if not isinstance(uniqSeedSeqs, (list, tuple)):
        raise ValueError("Input 'uniqSeedSeqs' must be a list or tuple.")
    if not all(isinstance(seq, str) for seq in uniqSeedSeqs):
        raise ValueError("Input 'uniqSeedSeqs' must contain only strings.")  
    lntFrequencies = []
    n = 0
    while n <= seedEnd - seedBegin:
        aCount = 0
        cCount = 0
        gCount = 0
        tCount = 0
        
        for uniqSeedSeq in uniqSeedSeqs:  
            nt = uniqSeedSeq[n]
            if   nt == 'a' or nt == 'A':
                aCount += 1
            elif nt == 'c' or nt == 'C':
                cCount += 1
            elif nt == 'g' or nt == 'G':
                gCount += 1
            elif nt == 't' or nt == 'T':
                tCount += 1                
            elif nt == 'u' or nt == 'U':
                tCount += 1    
                            
        n+=1
        ntCount = aCount + cCount + gCount + tCount
        lntFrequencies.append({'A': aCount/ntCount, 'C': cCount/ntCount, 'G': gCount/ntCount, 'T': tCount/ntCount})

    # convert list to dataframe
    dfNTFrequencies = pd.DataFrame(lntFrequencies)
    
    return dfNTFrequencies

def generateLogoPlot(dfNTFrequencies):
    import logomaker as lm
    import matplotlib.pyplot as plt
    
    logo = lm.Logo(dfNTFrequencies, font_name = 'Arial Rounded MT Bold')
    foldername = os.path.dirname(filename)
    basename = Path(filename).stem    
    outputpngfile = os.path.join(foldername, basename + "__uniqseeds_logoplt" + ".png")       
    plt.savefig(outputpngfile)

#
# def calcAverageGCPercent():
#     '''
#     calculate GC percent for each sequence and return the average value
#     :return:
#     '''
#     totalGCPercent = 0
#     sCount = 0
#     for seqLine in sequenceLines:
#         seq = sequence.Sequence(headerLines[sCount], seqLine)
#
#         seq.calcGC()
#
#         print("for sequence <" + seq.getHeaderLine() + "> GC% is <" + str(100.0*seq.getGCPercent()) + ">")
#         totalGCPercent = totalGCPercent + seq.getGCPercent()
#
#     return totalGCPercent/len(sequenceLines)
  
def initLogger(md5string):

    ''' setup log file based on project name'''
    projectBaseName = ""

    projectBaseName = Path(fastaFile).stem

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, projectBaseName + "__" + dt_string + "__" + md5string +".log")
    handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(level=logging.DEBUG)

    fileh = logging.FileHandler(logfileName, 'a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(formatter)

    log = logging.getLogger()  # root logger
    log.setLevel(logging.DEBUG)
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    log.addHandler(fileh)      # set the new handler
    log.addHandler(handler)
    logging.info("+" + "*"*78 + "+")
    logging.info("project log file is <" + logfileName + ">")
    logging.info("+" + "*"*78 + "+")
    logging.debug("debug mode is on")




# def parseArgs(argv):
#     '''parse out Command line options.'''
#
#     program_name = os.path.basename(sys.argv[0])
#     program_version = "v%s" % __version__
#     program_build_date = str(__updated__)
#     program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
#     #program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
#     program_license = '''%s
#     i
#       Created by Simon Rayner on %s.
#       Copyright 2023 Oslo University Hospital. All rights reserved.
#
#       Licensed under the Apache License 2.0
#       http://www.apache.org/licenses/LICENSE-2.0
#
#       Distributed on an "AS IS" basis without warranties
#       or conditions of any kind, either express or implied.
#
#     USAGE
#     ''' % (program_name, str(__date__))
#
#     try:
#         # Setup argument parser
#         parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
#         parser.add_argument("-f", "--fasta_file", dest="fastafile", action="store", help="fasta file for which you want to calc GC% [default: %(default)s]")
#
#         # Process arguments
#         args = parser.parse_args()
#
#         global fastaFile
#         global speciesCode
#         global seedBegin
#         global seedEnd
#
#         fastaFile = args.fastafile
#         speciesCode = args.speciescode
#         seedBegin = str(args.seedbegin)
#         seedEnd = str(args.seedend)
#
#
#         # check the user specified a fasta file, if not warn and and exit
#         if fastaFile:
#             print("fasta file is <" + fastaFile + ">")
#         else:
#             print("you must specify a fasta file")
#             exit
#
#         # Species Code is not required
#         if speciesCode:
#             print("speciesCode is <" + speciesCode + ">")
#
#         # check the user specified a start position for the seed region, if not warn and and exit
#         if seedBegin:
#             print("seed region start is <" + str(seedBegin) + ">")
#         else:
#             print("you must specify a start position for the seed region using the -b/--seed_begin flag")
#             exit
#
#         # check the user specified a stop position for the seed region, if not warn and and exit
#         if seedEnd:
#             print("seed region stop is <" + str(seedEnd) + ">")
#         else:
#             print("you must specify a stop position for the seed region using the -e/--seed_end flag")
#             exit
#
#
#     except KeyboardInterrupt:
#         ### handle keyboard interrupt ###
#         return 0
#     except Exception as e:
#         print(e)
#         if DEBUG or TESTRUN:
#             raise(e)
#         indent = len(program_name) * " "
#         sys.stderr.write(program_name + ": " + repr(e) + "\n")
#         sys.stderr.write(indent + "  for help use --help")
#         return 2


def calcAverageGCPercent():
    '''
    calculate GC percent for each sequence and return the average value
    :return:
    '''
    totalGCPercent = 0
    sCount = 0
    for seqLine in sequenceLines:
        seq = sequence.Sequence(headerLines[sCount], seqLine)

        seq.calcGC()

        print("for sequence <" + seq.getHeaderLine() + "> GC% is <" + str(100.0*seq.getGCPercent()) + ">")
        totalGCPercent = totalGCPercent + seq.getGCPercent()

    return totalGCPercent/len(sequenceLines)

    
    
    

def main(argv=None): # IGNORE:C0111

    if argv is None:
        argv = sys.argv

    #md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    #initLogger(md5String)
    speciesCode = "rno"
    global filename
    global seedBegin
    global seedEnd
    seedBegin = 2
    seedEnd = 8
    
    filename = r"C:\Users\user\AdvancedPrograming\day1\data\mature.fa"
    
    n = readFastaFile(filename)

    print (str(n))

    uniqSeedSeqs = getUniqueSeedSequences()

    print (uniqSeedSeqs)

    print('done')

    uniqSeedSeqs = writeUniqSeqs(uniqSeedSeqs)

    print(uniqSeedSeqs)
  

    # Call the function with the list of sequences

    dfNTFrequencies = getNucleotideFrequencyMatrix(uniqSeedSeqs)

    # Print the result

    print(dfNTFrequencies)

    generateLogoPlot(dfNTFrequencies)
     

if __name__ == '__main__':

    sys.exit(main())