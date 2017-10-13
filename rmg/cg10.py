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
                print("%12.3f " % (_time - timeref),end='')
                
        print(msg)
    
def line(level):
    if args.verbose >level:
        rows, columns = os.popen('stty size', 'r').read().split()
        print("-"*(int(columns)))
    

def build_models(args,inst,cg ):

    x0 = inst.map_assign()
    cg.build_model()

    cg.write_suffix("empty")
    
    for m in inst.M:
        msg(1,"building machine %5d (of %5d) MIP model" % (m,inst.nmach))
        cg.build_column_model(m)
            
        mres = inst.mach_validate(m,x0[:,m])
        cres = CG.CGColumn(obj=mres.obj,
                           procs = x0[:,m],
                           rc = 0,
                           g = None,
                           servs = None,
                           z = None,
                           hsigma = None,
                           ggamma = None,
                           pixp = None,
                           u = None,
                           ut = None,
                           a = None,
                           d = None,
                           b = None,
                           pmc = None,
                           mmc = None,
                           rtime = None
        )
                
        ares = cg.lp_add_col(m,cres)
        

def column_generation_loop(args,inst,cg):
    continue_cond = True

    if args.dump:
        cg.write_suffix("first_solve")
    
    res = cg.solve()
    first_obj = res.obj
    best_int_obj = first_obj

    stab = CG.Stabilization(inst, args)

    alpha = stab.alpha0()

    loop_start = time()
    while continue_cond:
        rtime = 0
        line(1)
        

        r_start = time()
        if args.dump:
            cg.write_suffix(stab.iterations())

        res = cg.solve()
        if args.dump:
            cg.writesol(stab.iterations())
        rtime+=res.rtime
        r_stop = time()

        msg(1,"%5d master %20.3f %20.3f %8.3f (%8.3f)   %.6f %20.3f %s" %(stab.iterations() ,res.obj,first_obj, r_stop - r_start, res.rtime, alpha, res.obj - first_obj, res.allint), timeref=loop_start)


        if res is None:
            raise Exception("LP não ótimo")
        
        rtime += res.rtime

        stduals = stab.alpha(res)
        omega = 0
        add_col = False
        for m in inst.M:
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
            rtime += cres.rtime
            
            ares = cg.lp_add_col(m,cres)
            
            if ares.status == CG.CGAddStatus.Added:
                add_col = True
            elif ares.status == CG.CGAddStatus.NotAdded:
                print(ares)
                raise Exception("problema ao adicionar coluna")

            c_stop = time()
            msg(3,"%5d  %5d %20.3f %20.3f %8.3f (%8.3f) %s" % ( stab.iterations(), m, cres.rc, cres.obj, c_stop -c_start ,cres.rtime ,ares.status.value), timeref=loop_start)
            
        if omega > - args.epslon or not add_col:
            continue_cond = False
        if stab.iterations() == args.maxsteps:
            continue_cond = False

        msg(1,"%5d  omega %20.3f %20.3f %8.3f (%8.3f)   %0.6f %5d %5d %8.3f %s" %(stab.iterations(),omega,stab.best_omega(), time() - r_start , rtime,alpha, stab.improvements(), stab.nonimprovements(),omega/stab.best_omega(), continue_cond),timeref=loop_start)
        
        alpha = stab.compute(omega = omega,stabdual = stduals)



def cuts_loop(args,inst,cg):
    line(1)
    _t = time()
    msg(1,"   Extending", timeref=_t)
    cg.extend()
    msg(1,"   Done, pre-gen all PPP cuts", timeref=_t)
    cg.cuts_prestats()
    msg(1,"   Done, first solve", timeref=_t)
    res = cg.solve()
    if res.allint:
        msg(1,"  LP solution is alreay integer, skipping all", timeref=_t) 
        return
    #cg.cuts_print_violated()
    
    #continue_cond = True
    #while continue_cond:
        
    



def mip(args,inst,cg):
    line(1)
    mip_start=time()
    msg(1,"   Converting ", timeref=mip_start)

    cg.convert()
    msg(1,"   Solving", timeref=mip_start)
    
    solution = cg.final_solve()

    line(1)
    msg(1,"Solution: %20.3f" % (solution.obj),timeref=mip_start)
    print(' '.join(str(i) for i in inst.assign()))
    print(' '.join(str(i) for i in solution.assign))

    if args.new_solution_filename:
        args.new_solution_filename.write(' '.join(str(i) for i in solution.assign))

    line(1)
    
def main():

    if args.name:
        print("Rodrigo Mosconi (1512344) <rmosconi@inf.puc-rio.br>")
        return

    global all_start
    all_start = time()
    msg(2,"Load Instance data")
            
    inst = Instance(args)

    msg(2,"CG Instance ")

    cg = CG.CG6(instance=inst,args=args)

    if args.loadfile:
        msg(0,"load saved model")
        cg.load(args.loadfile)
    else:
        msg(0,"building LP Model")

        build_models(args=args, inst=inst, cg=cg)

        column_generation_loop(args=args, inst=inst, cg=cg)
        if args.loadfile:
            msg(0,"saving the model")
            cg.save(args.savefile)


    cuts_loop(args=args, inst=inst, cg=cg)
    
    mip(args=args, inst=inst, cg=cg)



"""            
    if k >= 4 and False:
        continue_cond = False

    if continue_cond:
        if args.verbose > 1:
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
