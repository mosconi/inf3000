#!/usr/bin/env python3

import argparse,os,sys
import numpy as np

np.set_printoptions(formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+17d" % x, })

from time import time

from rmg import roadef,common,instance,Instance

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser
]
)

args = parser.parse_args()

if args.verbose is None:
    args.verbose = 1

all_start =0

def do_branches(instance, machine, assign, proc):
    nproc = instance.nproc
    P = instance.P
    nmach = instance.nmach
    M = instance.M

    for p in P:
        if p == proc: continue
        assign[p] = not assign[p]
        val = instance.mach_validate(machine = machine, map_assign = assign)
        if not 
            

def branch_and_bound(instance):
    nproc = instance.nproc
    P = instance.P
    nmach = instance.nmach
    M = instance.M

    # Adiciona mapa/lista de lista com todos os x validos?
    
    for m in M:
        x0 = instance.mach_map_assign(m).copy()
        x = x0
        for p in P:
            x[p] = not x[p]
            val = instance.mach_validate(machine = m, map_assign = x)
            if val.status:
                do_branches(instance = instance,
                            machine = m,
                            assign = x,
                            proc = p
                            )
                            
        
    

def main():
    
    if args.name:
        print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
        return

    global all_start
    all_start = time()
    msg(2,"Load Instance data")
            
    inst = Instance(args)
    branch_and_bound 



