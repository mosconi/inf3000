#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
from time import time

from gurobipy import *

def usage():
    print("%s [-n name] [-o outputfile ] [-s] [-q] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
name="roadef"
outputfile=None
savemodel=False
verbose=True
tex=True

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

if not os.path.exists(modelfile):
    sys.exit("File %s not found" % modelfile)

T=[]
Wlc=[]
M=[]
C=[]
C_bar=[]
S=[]
R=[]
RHO=[]
MU=[]

f = open(modelfile, "r")

lines = [list(map(int,l.rstrip('\n').split())) for l in f]

nres = lines.pop(0)[0]

for r in range(nres):
    l = lines.pop(0)
    T.append(l[0])
    Wlc.append(l[1])


Wlc=np.array(Wlc)
nmach = lines.pop(0)[0]

L=[[] for i in range(nmach)]
N=[[] for i in range(nmach)]

for m in range(nmach):
    l = lines.pop(0)
    N[l.pop(0)].append(m) # neighborhood
    L[l.pop(0)].append(m) # location
    C.append(l[:nres])
    del(l[:nres])
    C_bar.append(l[:nres])
    del(l[:nres])
    MU.append(l)

# remove empty neighborhood/locations
L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
C_bar=np.array(C_bar)
MU=np.array(MU)

nserv = lines.pop(0)[0]
delta=[0 for s in range(nserv)]
sdep={}

for s in range(nserv):
    l = lines.pop(0)
    delta[s]=l[0]
    if l[1]>0:
        sdep[s]=l[2:]

nproc = lines.pop(0)[0]

S=[[] for i in range(nserv)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    RHO.append(l[0])

nbal = lines.pop(0)[0]
bT=np.zeros((nres,nres),dtype=np.int32)
Wbal=np.zeros((nres,nres),dtype=np.int32)
for b in range(nbal):
    l = lines.pop(0)
    r1=l[0]
    r2=l[1]
    bT[ r1,r2]=l[2]
    l = lines.pop(0)
    Wbal[r1,r2] = l[0]

WPMC=lines[0][0]
WSMC=lines[0][1]
WMMC=lines[0][2]
    
R = np.array(R)
S = [s for s in S if s]    

x_bar= np.zeros((nproc,nmach), dtype=np.int32)


T=np.array(T)

if verbose: print(">> problem size")
if verbose: print(">>> resources ",nres)
if verbose: print(">>> transisent ",T.sum())
if verbose: print(">>> machines ",nmach)
if verbose: print(">>> process ",nproc)
if verbose: print(">>> services ",nserv)
if verbose: print(">>> balances ",nbal)
if verbose: print(">>> neighborhood ",len(N))
if verbose: print(">>> locations ",len(L))

if tex:
    print("\\verb|%s| & %d (%d/%d) & %d & %d  & %d & %d  & %d \\\\"  %
    (name,nres,T.sum(),nbal,nmach,nproc,nserv,len(N),len(L)))

    
