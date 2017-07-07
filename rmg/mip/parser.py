import argparse

parser = argparse.ArgumentParser(add_help=False,
    description="Arguments required for model/solver")

solver = parser.add_argument_group("solver","Options to change solver behaviour")
solver.add_argument("--epslon",type=float,default=10**-1,
                    help="Error tolerance")
solver.add_argument("--runname",dest="run_name",default="mip",
                    help="Name for the model")
solver.add_argument("--validate",dest="validate",
                    help="Validate solver results.")
log_group = parser.add_argument_group("log","For solver logging")
log_group.add_argument("--log",dest="log",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--console",dest="console",default=False,action="store_true",
                       help="Write solver messages to log")
log_group.add_argument("--logfile",dest="logfile",action="store",
                       help="File to Log")
