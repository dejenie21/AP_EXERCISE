'''
Created on Apr 22, 2023
@author: Dejenie Shiferaw
'''
#For the GC calc project: Add a logger method, generate logo plots and seed sequence from miRNAs for  mouse, human and rat by changing speciesCode.

import json
import sys
import os
import miRNA

def readFastaFile(filename):

    '''
    load specified fasta file and store header and sequence as entries in two lists
    :param self:
    :return:
    '''
    speciesCode = "hsa" #add species code
    print("load sequences from fasta file <" + filename+ ">")

    global headerLines
    global sequenceLines

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
  
    print("loaded <" + str(seqCount) + "> sequences and kept <"+ str(len(headerLines)) + "> with species code [" + speciesCode + "]")

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
        file.write(">uniqseed_" + str(s) + "\n")
        file.write(uSeq + "\n")

        s+=1

    file.close()    

    return uSeqs
    

def getNucleotideFrequencyMatrix(uniqSeedSeqs):
    '''
    calculate nucleotide frequencies at each position
    '''
    import numpy as np
    import pandas as pd    

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
            if nt == 'a' or nt == 'A':
                aCount += 1
            elif nt == 'c' or nt == 'C':
                cCount += 1
            elif nt == 'g' or nt == 'G':
                gCount += 1
            elif nt == 't' or nt == 'T':
                tCount += 1                
            elif nt == 'u' or nt == 'U':
                tCount += 1                             
        n += 1
        ntCount = aCount + cCount + gCount + tCount
        lntFrequencies.append({'A': aCount/ntCount, 'C': cCount/ntCount, 'G': gCount/ntCount, 'T': tCount/ntCount})

  # convert list to dataframe
    dfNTFrequencies = pd.DataFrame(lntFrequencies)
    return dfNTFrequencies

def generateLogoPlot(dfNTFrequencies):

    import logomaker as lm
    import matplotlib.pyplot as plt
    from pathlib import Path
    logo = lm.Logo(dfNTFrequencies, font_name = 'Arial Rounded MT Bold')
    foldername = os.path.dirname(filename)
    basename = Path(filename).stem    
    outputpngfile = os.path.join(foldername, basename + "__uniqseeds_logoplt" + ".png")       
    plt.savefig(outputpngfile)
    
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
   
def parseArgs():
    
    try:
        # Setup argument parser
        parser = ArgumentParser(description="parse test", formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-y", "--yaml_file", dest="yamlfile", action="store", help="yaml file name [default: %(default)s]")

        global yamlFile

        # Process arguments
        args = parser.parse_args()        
        yamlFile = args.yamlfile
        
    except Exception as e:
        print(e)
        return 2  

from pathlib import Path

from datetime import datetime
import hashlib
import logging   
    
def initLogger(md5string):

    ''' setup log file based on project name'''
    logfileBasename = ""

    logfileBasename = Path(yamlFile).stem

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, logfileBasename + "__" + dt_string + "__" + md5string +".log")
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


def main(argv=None):
    if argv is None:
        argv = sys.argv

    global speciesCode
    speciesCode = "hsa"
    global filename
    global seedBegin
    global seedEnd
    seedBegin = 2
    seedEnd = 8

    filename = r"C:\Users\user\AP\day6\data\five.fa"
    n = readFastaFile(filename)
    print (str(n))
    uniqSeedSeqs = getUniqueSeedSequences()
    print (uniqSeedSeqs)
    print('done')
    uniqSeedSeqs = writeUniqSeqs(uniqSeedSeqs)
    print(uniqSeedSeqs)
  
    dfNTFrequencies = getNucleotideFrequencyMatrix(uniqSeedSeqs)
    print(dfNTFrequencies)
    generateLogoPlot(dfNTFrequencies)
     

if __name__ == '__main__':

    sys.exit(main())
