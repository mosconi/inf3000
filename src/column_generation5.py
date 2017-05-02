#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from instance import Instance

from gurobipy import *

from columngeneration import CG5

rows, columns = os.popen('stty size', 'r').read().split()
np.set_printoptions(linewidth=int(columns)-5, formatter={'float_kind': lambda x: "%+20.3f" % x,'int': lambda x: "%+20d" % x, })

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=True
epslon = 10**-1
xi = 10**3
hsz = 10
slack = 0.02

try:
    opts, args = getopt.getopt(sys.argv[1:],"tqsn:m:a:o:h",["tex","quiet","save","name=","model=","assign=","output=","help"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit()
    
for o,a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ( "-a", "--assign"):
        assignfile = a
    elif o in ( "-m", "--model"):
        modelfile = a
    elif o in ( "-n", "--name"):
        name = a
    elif o in ( "-o", "--output"):
        outputfile = a
    elif o in ( "-s", "--save"):
        savemodel=True
    elif o in ( "-q", "--quiet"):
        verbose=False
    elif o in ( "-t", "--tex"):
        tex=True
    else:
        assert False, "unhandled option"
    
if modelfile is None:
    sys.exit("Model file not defined")

if assignfile is None:
    sys.exit("Assign file not defined")

if not os.path.exists(modelfile):
    sys.exit("File %s not found" % modelfile)

if not os.path.exists(assignfile):
    sys.exit("File %s not found" % assignfile)


inst = Instance(model=modelfile, assign=assignfile)


cg = CG5(inst)

#dir(cg)
cg.dual_history()
