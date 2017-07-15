import argparse
"""
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')

args = parser.parse_args()
print args.accumulate(args.integers)
"""
import sys
parser = argparse.ArgumentParser(description='###')
parser.add_argument('ofile')
parser.add_argument('--n',type=int,default=6,dest='n')
#parser.add_argument('--td',type=int, default=0,dest='timedelay')
args = parser.parse_args()

#import time
#time.sleep(args.timedelay)

import numpy as np
a = np.arange(args.n)*np.ones((7,args.n))

np.savez(args.ofile,a)