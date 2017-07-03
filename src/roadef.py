import argparse

import numpy as np
    
parser = argparse.ArgumentParser(prog="ROADEF",add_help=False,
                                 description="Arguments required by roadef")

# Descritos no roadef
roadef = parser.add_argument_group("roadef","Requested arguments by the challenge")
roadef.add_argument("-name",
                    action="store_true",
                    help="Return the identifier of the team that is the author of the executable. *In this case, the author name*")

roadef.add_argument("-t",dest="timelimit", type=int, default=np.inf,
                    help="Stop the program execution after time_limit seconds")
roadef.add_argument("-p",dest="instance_filename", type=argparse.FileType('r'),
                    help="Load the data associated with the instance instance_filename")
roadef.add_argument("-i",dest="original_solution_filename", type=argparse.FileType('r'),
                    help="Designate the file with the reference solution")
roadef.add_argument("-o",dest="new_solution_filename", type=argparse.FileType('wb',0),
                    help="Designate the result file")
roadef.add_argument("-s",dest="seed",
                    help="Force program with random to be deterministic")

args = parser.parse_known_args()[0]
