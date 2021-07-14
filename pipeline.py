#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging


# Third party imports
import argparse
import subprocess
import datetime
import glob2


# Local application imports
from . import version


"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com), Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 06 July 2021

TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

# COLORS AND AND FORMATTING

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

logger = logging.getLogger()


def run_submodule(parser, args):
    if args.command == 'basecall':
        from . import guppy_minion as submodule
    if args.command == 'variants':
        from . import variants as submodule
    if args.command == 'assembly':
        from . import assembly as submodule

    # Run the chosen submodule
    submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument('-q', '--quiet', help = 'Do not output warnings to stderr', action = 'store_true', dest = 'quiet')


def init_pipeline_parser():
    """Wraps the argparse parser initialisation
    Returns
    -------
    argparse.ArgumentParser
        The initialised argparse Argument Parser for the pipeline
    """

    parser = argparse.ArgumentParser(prog = 'prokaion', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', help = 'Installed Prokaion version', action = 'version', version = '%(prog)s' + str(version.__version__))

    subparsers = parser.add_subparsers(title = '[Sub-Commands]', dest = 'Command', parser_class = ArgumentParserWithDefaults)



def main():

    # Init the pipeline parser
    parser = init_pipeline_parser()

    # Collect the args from submodules
    args = parser.parse_args(sys.argv[1:])

    if args.quiet:
        logger.setLevel(logging.ERROR)

    # Run the subcommand or print usage if no subcommand provided
    if args.command:
        args.func(parser, args)
    else:
        parser.print_usage()



if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise