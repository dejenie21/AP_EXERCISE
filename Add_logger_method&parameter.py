'''
Created on Apr 22, 2023
@author: Dejenie Shiferaw
'''
#Add a logger method to the complexity project
#Add a parameter using parseargs that sets nmax

from timeit import default_timer as timer 
from argparse import ArgumentParser 
from argparse import RawDescriptionHelpFormatter 
from datetime import datetime

import unittest 
import sys 
import os 
import hashlib
from pathlib import Path

import logging
logger = logging.getLogger(__name__) 


def initLogger(md5string):

    now = datetime.now()
    dt_string = now.strftime("%Y%m%d_%H%M%S")
    logFolder = os.path.join(os.getcwd(), "logfiles")
    if not os.path.exists(logFolder):
        print("--log folder <" + logFolder + "> doesn't exist, creating")
        os.makedirs(logFolder)
    logfileName = os.path.join(logFolder, "__" + dt_string + "__" + md5string +".log")
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


# def parseArgs():
#     '''parse out Command line options.'''
#     try:
#         parser = ArgumentParser(description="a program to calculate Fibonacci's number", formatter_class=RawDescriptionHelpFormatter)
#         parser.add_argument("-m", "--max_number", dest="maxnumber", action="store", help="max value to calculate Fibonacci number [default: %(default)s]")
#         parser.add_argument("-n", "--nmax", dest="nmax", action="store", help="max number of Fibonacci to be calculated [default: %(default)s]")
#
#         # Process arguments
#         args = parser.parse_args()
#         global maxNumber, nmax
#         maxNumber = args.maxnumber
#         nmax = int(args.nmax) if args.nmax else 10
#
#      # Set maxNumber to 10 if not specified
#         if not maxNumber:
#             maxNumber = 10
#         # check the user specified a max number of Fibonacci, if not warn and and exit
#         if maxNumber:
#             logger.info("max Fibonacci number to calculate is <" + str(maxNumber) + ">")
#         else:
#             logger.error("you must specify a max Fibonacci number")
#             sys.exit(1)
#     except KeyboardInterrupt:
#         ### handle keyboard interrupt ###
#         return 0
#     except Exception as e:
#         logger.error(e) 
def timeSimpleVersusDynamicFibo(nMax):
    '''
    compare times for simple and dynamic methods to calculate Fibonacci's number
    return fibonacciTimes
    '''
    n=0
    fibonacciTimes=[]
    while n<nMax:
        logger.debug("--- " + str(n))
        simpleStartTime = timer()
        fibSimple(n)
        simpleEndTime = timer()
        dynamicStartTime = timer()
        fibDynamic(n)
        dynamicEndTime = timer()
        fibonacciTimes.append({"n":n, "simple": (simpleEndTime-simpleStartTime), "dynamic": (dynamicEndTime-dynamicStartTime) })
        n+=1
    return fibonacciTimes
def generateTimingPlot(fibonacciTimes, nmax):
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    dfFibonacciTimes = pd.DataFrame(fibonacciTimes)
    dfMelt = dfFibonacciTimes.melt(id_vars=['n'], value_vars=['simple','dynamic'])
    #sns.lineplot(data=dfMelt, x="n", y="simple", hue='variable')
    plt.xlabel("Fibonacci number")
    plt.ylabel("runtime (s)",)
    sns.scatterplot(data=dfMelt, x="n", y="value", hue='variable').set(title='Title of Plot')
    timingPlotFile = os.path.join(os.getcwd(), "fibonacci_timing_0_to_" + str(nmax) + ".png")
    plt.savefig(timingPlotFile) 
def fibSimple(n):
    # base case
    if n == 0:
        return(0)
    if n == 1:
        return(1)
    # pattern: ith = (i-1)th + (i-2)th
    return(fibSimple(n-1) + fibSimple(n-2)) 
def writeTimingDataToFile(fiboTimes, nmax):  
    import csv

    timingDataFile = os.path.join(os.getcwd(), "fibonacci_timing_0_to_" + str(nmax) + ".tsv")
    print("timing data will be written to <" + timingDataFile + ">")
    with open(timingDataFile, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(['n', 'recursive', 'dynamic'])
        for timing in fiboTimes:
            writer.writerow(timing.values())    
    
    print("---done")
    
def fibDynamic(n):

    
    # base case
    if n == 0:
        return(0)
    if n == 1:
        savedFibNumbers[1] = 1
        return(1)
    
    if savedFibNumbers[n] != 0:
        return savedFibNumbers[n]
    
    savedFibNumbers[n] = fibDynamic(n -1) + fibDynamic(n - 2)

    return(savedFibNumbers[n])    
 
class CheckFibonacciDynamic(unittest.TestCase):
    
    # check the fibonacci dynamic calculation is correct
    # can read more about unittesting here 
    #  https://machinelearningmastery.com/a-gentle-introduction-to-unit-testing-in-python/
    def test_negative(self):
        
        nmax= 10
        global savedFibNumbers
        savedFibNumbers=[0]*nmax
        
        message = "fibonacci calculation is wrong !"
        
        # assertEqual() to check equality of first & second value

        self.assertEqual(fibDynamic(10), 55, message)
  
        

def main(argv=None): 

    if argv is None:
        argv = sys.argv
    nmax = 10
    
    # parseArgs()
    
    global savedFibNumbers
    savedFibNumbers=[0]*nmax
    

    #unittest.CheckFibonacciDynamic()
    
    fibonacciTimes = timeSimpleVersusDynamicFibo(nmax)
    writeTimingDataToFile(fibonacciTimes, nmax)
    generateTimingPlot(fibonacciTimes, nmax)
           
    md5String = hashlib.md5(b"CBGAMGOUS").hexdigest()
    # parseArgs()
    
    initLogger(md5String)
    
if __name__ == '__main__':

    sys.exit(main())

