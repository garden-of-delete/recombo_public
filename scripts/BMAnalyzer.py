import sys
import numpy as np
import random
import argparse
#set up argument parser
parser = argparse.ArgumentParser(description='Analyzes and blocks a sequence of values drawn from a MCMC process into independent batches.')
parser.add_argument('infile',nargs='+',help='sequence of values separated by a newline character',)
parser.add_argument('outfile',nargs='+',help='name of the output file',)
parser.add_argument('--start',help='specify the starting block size. 1 is default and assumes independence.',dest='start',type=int)
parser.add_argument('--isize',help='specifies the amount to increase the block size by with each iteration. 1 is deafult.',dest='isize',type=int)
parser.set_defaults(start=1,isize=1)
args = parser.parse_args()

class BMAnalyzer(object):
    #constructor
    def __init__(self,chain_data):
        #elements for use in block mean analysis
        self.chain_data = chain_data #stores the input data
        self.n_lines = len(chain_data) #number of lines in the .info file
        self.block_size = 0
        self.n_blocks = 0 #integer floor of n_lines/block_size

        return super(BMAnalyzer, self).__init__()

    def write_datafile(self,filename):
        datafile = open(filename,'w')
        for item in self.chain_data:
            datafile.write(repr(item)+'\n')
        datafile.close()

    def compute_block_means(self,block_size,input_list):
        """takes a n input list and a block size, forming a new list of lists, each sublist representing one block"""
        # group the list into blocks of size block_size
        if block_size == 1:
            return np.asarray(input_list, dtype=np.float64)
        block_index = 0
        sum = 0
        block_means = []
        #for value in input_list:
        for i in range(len(input_list)):
            if block_index < block_size:
                sum += input_list[i]#value
                block_index += 1
            else:
                sum = np.float64(sum) / np.float64(block_size)  # 64bit floating point precision
                block_means.append(sum)  # 32bit floating point precision
                # block_index = 0
                sum = 0
                sum += input_list[i]#value
                block_index = 1
        if block_index == block_size:
            sum = np.float64(sum) / np.float64(block_size)  # 64bit floating point precision
            block_means.append(sum)  # 32bit floating point precision
        return np.asarray(block_means, dtype=np.float64)

    def block_length_analysis(self,starting_block_size,iter_size,debug):
        """perform block-mean analysis on the length of the chain with the largest z-value"""
        # perform block mean analysis. stores block size and associated results to index.block_size.
        #print "Performing block mean analysis on " + data.name + "..."
        # start with block_size 1 (assumption of independent data)
        self.block_size = starting_block_size
        test = False
        # while statistical test has not passed yet
        while test == False:
            # calculate block means
            if debug:
                print 'Current blocksize: ' + repr(self.block_size)
                print 'Chain data length: ' + repr(len(self.chain_data))
            # block_means = self.compute_block_means(topo.block_size, [element[0] for element in topo.chain_data[0]])
            block_means = self.compute_block_means(self.block_size, self.chain_data)
            if debug:
                print'Blocked list size: ' + repr(len(block_means))
            self.n_blocks = self.n_lines / self.block_size # number of full blocks. partial block on the end is discarded.
            #ensure we have enough blocks for statistical significance
            if self.n_blocks < 30:
                print 'Error: Not enough blocks. Provide more data.'
                return
            # calculate the numerator of the test statistic
            num = 0.0
            for i in range(self.n_blocks - 1):
                num += (block_means[i] - block_means[i + 1]) ** 2 #possible bug here
            # calculate variance of block means
            var = np.var(block_means)
            # calculate denominator of test statistic
            den = 2 * (self.n_blocks - 1) * (var)
            print num,den
            # calculate test statistic
            if debug:
                print float((self.n_blocks - 2)), float((self.n_blocks ** 2 - 1)), float((self.n_blocks - 2)) / float(
                    (self.n_blocks ** 2 - 1)), np.sqrt(float((self.n_blocks - 2)) / float((self.n_blocks ** 2 - 1)))
            test_stat = 1 - float(num) / float(den)
            comparison = 1.645 * np.sqrt(float((self.n_blocks - 2)) / float((self.n_blocks ** 2 - 1)))
            # if in debug mode, print the test statistic and comparison value
            if debug:
                print 'Null Hypothesis: ' + repr(test_stat)
                print 'Comparison: ' + repr(comparison)
                print ''
            # perform test
            if test_stat <= comparison:
                test = True
            else:
                self.block_size += iter_size
        print "... Complete. Block size is: " + repr(self.block_size)

#read datafile
def read_datafile(filename):
    temp = []
    datafile = open(filename,'r')
    for entry in datafile:
        temp.append(float(entry.strip()))
    temp
    datafile.close()
    return temp

#DEMO: generate sequence of values from random 1d integer walk
def oneDwalk(initial_value,step_size,n_samples,seed):
    walk = []
    current_value = initial_value
    random.seed(seed)
    for i in range(0,n_samples):
        for j in range(0,step_size):
            current_value += np.random.uniform(-1,1,1).tolist()[0]
        walk.append(current_value)
    return walk

#Main
data = read_datafile(sys.argv[1])
bma = BMAnalyzer(data)
bma.block_length_analysis(args.start,args.isize,True)
bma.write_datafile(sys.argv[2])