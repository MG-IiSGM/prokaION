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
from . import pipeline


if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main())