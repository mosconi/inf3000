import argparse

parser = argparse.ArgumentParser(add_help=False,
    description="Arguments required for model/solver")

solver = parser.add_argument_group("solver","Options to change solver behaviour")
solver.add_argument("--epslon",type=float,default=10**-1,
                    help="Error tolerance (column generation)")
solver.add_argument("--tol",type=float,default=10**-6,
                    help="Error tolerance (calculations)")
solver.add_argument("--runname",dest="run_name",default="columngeneration",
                    help="Name for the model")
solver.add_argument("--novalidate",dest="validate",action="store_false",default=True,
                    help="Validate solver results.")
solver.add_argument("--nodump",dest="dump",action="store_false",default=True,
                    help="Dump models.")
log_group = parser.add_argument_group("log","For solver logging")
log_group.add_argument("--log",dest="log",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--console",dest="console",default=False,action="store_true",
                       help="Write solver messages to console")
log_group.add_argument("--logfile",dest="logfile",action="store",default="columngeneration.log",
                       help="File to Log")
