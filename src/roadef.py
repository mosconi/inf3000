#!/usr/bin/env python3

import os, sys, getopt
import numpy as np
import pulp
import scipy.optimize as opt

def usage():
    print("%s -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None


try:
    opts, args = getopt.getopt(sys.argv[1:],"m:a:h",["model=","assign=","help"])
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

T=[]
Wlc=[]
M=[]
C=[]
SC=[]
S=[]
R=[]
PMC=[]
MMC=[]

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
    SC.append(l[:nres])
    del(l[:nres])
    M.append(l)

L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
SC=np.array(SC)
M=np.array(M)

nserv = lines.pop(0)[0]

for s in range(nres):
    l = lines.pop(0)

nproc = lines.pop(0)[0]
S=[[] for i in range(nproc)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    PMC.append(l[0])



S = [s for s in S if s]    

R=np.array(R)

x_bar= np.zeros((nproc,nmach), dtype=np.int32)
x= np.zeros((nproc,nmach), dtype=np.int32)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]

for p in range(len(line)):
    x_bar[p,line[p]] = 1

T=np.array(T)

u_bar=T.reshape((2,1))*R.T.dot(x_bar)

fun = lambda xi: np.sum(xi*Wlc)


