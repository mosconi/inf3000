import argparse

parser = argparse.ArgumentParser(add_help=False,
                                 description="Common arguments for all modules")


# Other parameters:
verbosity = parser.add_argument_group("verbosity","Debug messages level")

verbosity.add_argument("-q","--quiet",action="store_const",dest="verbose",const=0,
                       help="Don`t print messages")
verbosity.add_argument("-v","--verbose",action="count",
                       help="Increase ouput verbosity")


args= parser.parse_known_args()[0]
