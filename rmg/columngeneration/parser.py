import argparse

parser = argparse.ArgumentParser(add_help=False,
    description="Arguments required for model/solver")

solver = parser.add_argument_group("solver","Options to change solver behaviour")
solver.add_argument("--epslon",type=float,default=1e-1,
                    help="Error tolerance (column generation)")
solver.add_argument("--tol",type=float,default=1e-6,
                    help="Error tolerance (calculations)")
solver.add_argument("--validate",dest="validate",action="store_true",default=False,
                    help="Validate solver results.")
solver.add_argument("--accept",dest="accept",action="store_true",default=False,
                    help="accept solver results.")


log_group = parser.add_argument_group("log","For solver logging")
log_group.add_argument("--log",dest="log",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--console",dest="console",default=False,action="store_true",
                       help="Write solver messages to console")
log_group.add_argument("--logfile",dest="logfile",action="store",default="columngeneration.log",
                       help="File to Log")

model = parser.add_argument_group("model","Options to change model behaviour")
model.add_argument("--pregenerate",dest="generate",default=False,action="store_true",
                       help="Pre-generate some columns")
model.add_argument("--mipdump",dest="mipdump",action="store_true",default=False,
                    help="Dump MIP model.")
model.add_argument("--runname",dest="run_name",default="columngeneration",
                    help="Name for the model")
model.add_argument("--mipstats",dest="mipstats",default=False,action="store_true",
                    help="Print the MIP model stats")


stab = parser.add_argument_group("stabilization","Options to change dual stabilization")
stab.add_argument("--start",dest="alpha0",default=0.9,type=float,
                  help="Starting alpha")
stab.add_argument("--min",dest="alpha_min",default=0.0,type=float,
                  help="Ending alpha")
stab.add_argument("--method",dest="method",choices=["linearkp1","linearip1","exp","omegaratio","lineark","lineari","free","constant"],default="free",
                  help="Computation method")
stab.add_argument("--offset",dest="alpha_offset",default=0,type=int,
                  help="Offset to compute")
stab.add_argument("--scale",dest="alpha_scale",default=1.0,type=float,
                  help="Scale to compute")
stab.add_argument("--steps",dest="alpha_steps",default=1000,type=int,
                  help="Steps to compute")
stab.add_argument("--maxsteps",dest="maxsteps",default=-1,type=int,
                  help="MaxSteps to compute")
