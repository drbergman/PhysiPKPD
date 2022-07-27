#!/usr/bin/python

import sys, getopt, csv
import matplotlib.pyplot as plt
import numpy as np

def createPlots(inputfile, outputfile):
      
    x = []
    y = []
      
    try:
        with open(inputfile,'r') as csvfile:
            lines = csv.reader(csvfile, delimiter=',')
            for row in lines:
                x.append(row[0])
                y.append(int(row[1]))

        plt.plot(x, y)
        ax = plt.gca()
        ax.set_xticks(ax.get_xticks()[::100])
        #ax.set_xticks(np.arange(0, max(x), 10)) #step 10 digits
        plt.xlabel('Time (in minutes)')
        plt.ylabel('Cell Count')
        plt.title('Tumor Cell Count', fontsize = 20)
        #plt.show()
        plt.savefig(outputfile)
    except OSError as err:
        print("OS error: {0}".format(err))
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('createDataPlots.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('createDataPlots.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    createPlots(inputfile, outputfile)
    
if __name__ == "__main__":
    main(sys.argv[1:])