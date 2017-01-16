#!/usr/bin/env python

import os, sys, getopt
import numpy as np
from time import time

def usage():
    print("%s [-n name] -m modelfile -a assignfile\n" % sys.argv[0])

modelfile=None
assignfile=None
name="roadef"

try:
    opts, args = getopt.getopt(sys.argv[1:],"n:m:a:h",["name=","model=","assign=","help"])
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

# remove empty neighborhood/locations
L = [l for l in L if l]
N = [n for n in N if n]

    
C=np.array(C)
SC=np.array(SC)
MMC=np.array(M)

nserv = lines.pop(0)[0]

for s in range(nserv):
    l = lines.pop(0)

nproc = lines.pop(0)[0]

S=[[] for i in range(nserv)]

for p in range(nproc):
    l = lines.pop(0)
    S[l.pop(0)].append(p)
    R.append(l[:nres])
    del(l[:nres])
    PMC.append(l[0])

nbal = lines.pop(0)[0]
for b in range(nbal):
    l = lines.pop(0)
    l = lines.pop(0)
    
WPMC=lines[0][0]
WSMC=lines[0][1]
WMMC=lines[0][2]
    
R = np.array(R)
S = [s for s in S if s]    

x_bar= np.zeros((nproc,nmach), dtype=np.int32)

f = open(assignfile, "r")
line = [list(map(int,l.rstrip('\n').split())) for l in f][0]
assign=line[:]
for p in range(nproc):
    x_bar[p,line[p]] = 1

T=np.array(T)

print(">> problem size")
print(">>> resources ",nres)
print(">>> machines ",nmach)
print(">>> process ",nproc)
print(">>> services ",nserv)
print(">>> neighborhood ",len(N))
print(">>> locations ",len(L))
