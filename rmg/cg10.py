#!/usr/bin/env python3

import argparse,os,sys
import numpy as np

np.set_printoptions(formatter={'float_kind': lambda x: "%+17.3f" % x,'int': lambda x: "%+17d" % x, })

from time import time

from rmg import roadef,common,instance,Instance,columngeneration

import rmg.columngeneration as CG

parser = argparse.ArgumentParser(description="",
                                 parents=[roadef.parser,
                                          common.parser,
                                          instance.parser,
                                          columngeneration.parser
                                          
]
)


# Other parameters:
parser.add_argument("--save",dest="savefile",
                    help="Save a LP model")
parser.add_argument("--load",dest="loadfile",
                    help="Load a saved model")
args = parser.parse_args()

if args.verbose is None:
    args.verbose = 1


all_start =0
def msg(level, msg, timeref=None):
    global all_start 

    if args.verbose >level:
        if args.time:
            _time = time()
            print("%12.3f " % (_time  - all_start),end='')
            if timeref is not None:
                print("%12.3f " % (_time - t),end='')
                
        print(msg)
    
def line():
    rows, columns = os.popen('stty size', 'r').read().split()
    print("-"*(int(columns)-2))
    
        
def main():
    print("teste")
    if args.name:
        print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
        return

    global all_start
    all_start = time()
    msg(2,"Load Instance data")
            
    inst = Instance(args)

    cg = CG.CG6(instance=inst,args=args)

    if args.loadfile:
        msg(0,"load saved model")
        cg.load(args.loadfile)
    else:
        msg(0,"building LP Model")

    cg.build_model()

    for m in inst.M:
        msg(1,"building machine %5d (of %5d) MIP model" % (m,inst.nmach))
        cg.build_column_model(m)
        
    return



"""            
continue_cond = True

res = cg.solve_relax()
first_obj = res.obj
best_int_obj = first_obj

stab = CG.Stabilization(inst, args)

alpha = stab.alpha0()

k = 0
loop_start = time()
while continue_cond:
    rtime = 0
    _time = time()
    if args.verbose >1:
        print("-"*(int(columns)-2))
        if args.time:
            print("%12.3f %12.3f " % (_time - all_start, _time - loop_start),end='')
        print("%5d master" % (k),end=' ')
    if args.log:
        cg.lplog("\n/*"+"*"*70)
        cg.lplog(" *")
        cg.lplog(" * %12.3f %12.3f MASTER iteration %d" % (_time-all_start,_time-all_start,k))
        cg.lplog(" *")
        

    r_start = time()
    if args.dump:
        cg.lpwrite(k)

    res = cg.solve_relax()
    if args.dump:
        cg.lpwritesol(k = k)
    r_stop = time()

    if args.verbose >1:
        print("%20.3f %20.3f %8.3f (%8.3f)   %.6f %17.3f %s" % (res.obj,first_obj, r_stop - r_start, res.rtime, alpha, res.obj - first_obj, res.allint))
    if args.log:
        cg.lplog("\n/*"+"*"*70)
        cg.lplog(" *")
        cg.lplog(" *T:[MASTER:%05d] %12.3f %12.3f %20.3f %20.3f %8.3f (%8.3f) %20.3f" % (k,_time - all_start, _time - loop_start ,res.obj,first_obj, r_stop - r_start, res.rtime, res.obj - first_obj))
        cg.lplog(" *")

    #if  abs(first_obj - res.obj) > args.epslon: break

    rtime += res.rtime
    if res is None:
        raise Exception("LP não ótimo")
    

    stduals = stab.alpha(res)
    omega = 0
    add_col = False

    for m in range(inst.nmach):
        _time = time()
        if args.verbose>3:
            if args.time:
                print("%12.3f %12.3f " % (_time - all_start, _time - loop_start),end='')
            print("%5d  %5d" % (k,m),end=' ')
        if args.log:
            cg.machlog(m,"\n/*"+"*"*70)
            cg.machlog(m," *")
            cg.machlog(m," * %12.3f %12.3f machine[%d] iteration %d" % (_time-all_start,_time-all_start,m,k))
            cg.machlog(m," *")

        c_start = time()
        sigma = stduals.eta_lb[:,m] + stduals.eta_ub[:,inst.iN[m]] + stduals.omikron_lb[:,m] + stduals.omikron_ub[:,inst.iL[m]]

        cg.column_prepare(machine = m,
                          pi = stduals.pi,
                          gamma = stduals.gamma,
                          sigma = sigma,
                          mu = stduals.mu[m])
        if args.dump:
            cg.column_write(machine = m, k = k)
        
        cres = cg.column_compute(machine = m)

        c_stop = time()
        if cres is None:
            raise Exception("problema ao calcular coluna")


        if args.dump:
            cg.column_writesol(machine = m, k = k)

        if args.validate:
            vres = cg.column_validate(cres,
                                      machine = m,
                                      pi = stduals.pi,
                                      gamma = stduals.gamma,
                                      sigma = sigma,
                                      mu = stduals.mu[m])

            if vres.status != CG.CGValidateStatus.Valid:
                print(vres)
                raise Exception("problema ao validar coluna")
        
        omega += cres.rc
        ares = cg.lp_add_col(m,cres)

        if ares.status == CG.CGAddStatus.Added:
            add_col = True
            
        elif ares.status == CG.CGAddStatus.NotAdded:
            print(ares)
            raise Exception("problema ao adicionar coluna")
        
        if args.verbose>3:
            print("%20.3f %20.3f %8.3f (%8.3f) %s" %(cres.rc, cres.obj, c_stop -c_start ,cres.rtime ,ares.status.value),end=' ' )
            if ares.status == CG.CGAddStatus.Exist:
                print(ares.var.VarName)
            else:
                print()
        rtime += cres.rtime
        if args.log:
            cg.machlog(m,"\n/*"+"*"*70)
            cg.machlog(m," *")
            cg.machlog(m," *T:[%06d:%05d] %12.3f %12.3f %20.3f %20.3f %8.3f (%8.3f) %s" % (m,k,_time - all_start, _time - loop_start, cres.rc, cres.obj, c_stop -c_start ,cres.rtime ,ares.status.value))
            cg.machlog(m," *")

    _time= time()
                
    if args.verbose>1:
        if args.time:
            print("%12.3f %12.3f " % (_time - all_start, _time - loop_start),end='')
        print("%5d  omega %20.3f %20.3f %8.3f (%8.3f)" %(k,omega,stab.best_omega(), _time - r_start , rtime),end=' ')

    
    if omega > - args.epslon or not add_col:
        continue_cond = False
    if k >= 4 and False:
        continue_cond = False

    if continue_cond:
        if args.verbose > 1:
            print("C %0.6f %5d %5d %0.3f" %(alpha, stab.improvements(), stab.nonimprovements(),omega/stab.best_omega()))
        _c = "C"
    else:
        if args.verbose > 1:
            print("T %0.6f %5d %5d %0.3f" %(alpha, stab.improvements(), stab.nonimprovements(),omega/stab.best_omega()))
        _c = "T"
    if args.log:
        cg.lplog("\n/*"+"*"*70)
        cg.lplog(" *")
        cg.lplog(" *T:[ OMEGA:%05d] %12.3f %12.3f %20.3f %20.3f %8.3f (%8.3f) %s %.6f %d %d" % (k,_time - all_start, _time - loop_start, omega,stab.best_omega(), _time - r_start , rtime, _c,alpha, stab.improvements(), stab.nonimprovements()))
        cg.lplog(" *")


    alpha = stab.compute(omega = omega,stabdual = stduals)
    
    k+=1

    if (args.maxsteps > 0) and \
       (k > args.maxsteps):
        break
    if time() - loop_start > args.timelimit:
        break

    if res.obj < first_obj and res.obj - int(res.obj) < 1.0e-6:
        best_int_obj = int(res.obj)

#input('Press Enter to continue')


_time1 = time()
if args.verbose>1:
    print("-"*(int(columns)-2))
if args.verbose>1:
    if args.time:
        print("%12.3f " % (time() - all_start), end='')
    print("Extend - solve")

_time2 = time()
cg.extend()
res = cg.solve_relax()
_time3 = time()

if args.verbose >1:
    if args.time:
        print("%12.3f %12.3f" % (time() - all_start, _time3-_time1), end=' ')
    print("%20.3f %20.3f %8.3f (%8.3f %8.3f)   %.6f %17.3f" % (res.obj,first_obj, _time3 - _time1, _time3 - _time2, res.rtime, alpha, res.obj - first_obj))
        

if args.verbose>1:
    print("-"*(int(columns)-2))
input('Press Enter to continue')
if args.verbose>1:
    print("-"*(int(columns)-2))

cont_cond = True
c =0 
while cont_cond:
    c+=1
    if args.verbose>1:
        print("-"*(int(columns)-2))
        if args.time:
            print("%12.3f " % (time() - all_start), end='')
        print("solve - add")

    _time1 = time()
    res = cg.solve_relax()
    _time2 = time()
    if args.verbose >1:
        if args.time:
            print("%12.3f %12.3f" % (time() - all_start, _time2-_time1), end=' ')
        print("%20.3f %20.3f %8.3f (%8.3f %8.3f)   %.6f %17.3f %s" % (res.obj,first_obj, _time3 - _time1, _time3 - _time2, res.rtime, alpha, res.obj - first_obj, res.allint))

    try:
        while cg.cuts_add() is None:
            pass
        else:
            if args.verbose>1:
                print("-"*(int(columns)-2))
                if args.time:
                    print("%12.3f " % (time() - all_start), end='')
                    
                print("solve - filter")
                
            res = cg.solve_relax()
            _time3 = time()

        if args.verbose >1:
            if args.time:
                print("%12.3f %12.3f" % (time() - all_start, _time3-_time1), end=' ')
                print("%20.3f %20.3f %8.3f (%8.3f %8.3f)   %.6f %17.3f %s" % (res.obj,first_obj, _time3 - _time1, _time3 - _time2, res.rtime, alpha, res.obj - first_obj, res.allint))
        
        changed = cg.cuts_filter()
    except:
        break

    #if changed == 0:
    #    cont_cond = False
    if c >= inst.nproc*inst.nmach:
        break
    
if res.allint:
    cg.print_solution_relax()


if args.verbose>1:
    print("-"*(int(columns)-2))
    if args.time:
        print("%12.3f " % (time() - all_start), end='')
    print("Converting LP model to MIP model")
cg.lp2mip()
if args.verbose>1:
    if args.time:
        print("%12.3f " % (time() - all_start), end='')
    print("Done")

if args.mipstats:
    if args.verbose>1:
        print("-"*(int(columns)-2))
    if args.verbose>2:
        if args.time:
            print("%12.3f " % (time() - all_start), end='')
        print("MIP STATS")
    cg.mip_stats()

if args.dump or args.mipdump:
    cg.write()
    

mip_start=time()
#solution = cg.solve()

solution = cg.solve()

if args.verbose>0:
    if args.verbose>1:
        print("-"*(int(columns)-2))
    if args.time:
        print("%12.3f %12.3f Solution: %20.3f" % (time() - all_start,time() - mip_start,solution.obj))
    print(' '.join(str(i) for i in inst.assign()))
    print(' '.join(str(i) for i in solution.assign))
    if args.verbose>1:
        print("-"*(int(columns)-2))

if args.new_solution_filename:
    args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))
    
    if args.verbose>1:
        print("-"*(int(columns)-2))
    if args.time:
        _time = time()
        print("%12.3f Done in %12.3f/%12.3f s" % (_time - all_start,_time - loop_start,_time - mip_start))
"""
