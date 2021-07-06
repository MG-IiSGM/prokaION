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

from misc_ion import check_create_dir, check_file_exists, check_remove_file, execute_subprocess, check_reanalysis


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


def main():
    
    """
    Create main function to capture code errors: https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python
    """

    # ARGUMENTS

    def get_arguments():