import argparse

parser = argparse.ArgumentParser(add_help=False,
                                 description="Common arguments for all modules")


parser.add_argument("--time",dest="time",default=False,action="store_true",
                    help="Time the model")
parser.add_argument("--dump",dest="dump",action="store_true",default=False,
                    help="Dump models.")


# Other parameters:
verbosity = parser.add_argument_group("verbosity","Debug messages level")

verbosity.add_argument("-q","--quiet",action="store_const",dest="verbose",const=0,
                       help="Don`t print messages")
verbosity.add_argument("-v","--verbose",action="count",dest="verbose",default=1,
                       help="Increase ouput verbosity")


