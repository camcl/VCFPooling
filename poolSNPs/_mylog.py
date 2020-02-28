"""
Automatically outputs console output to a log-file
"""

import logging
import logging.handlers
import os
import sys

# handler = logging.handlers.WatchedFileHandler(
#     os.environ.get("LOGFILE", os.path.join(os.getcwd(),
#                                            os.path.basename(__file__).rstrip('.py'),
#                                            '.log')))
# formatter = logging.Formatter(logging.BASIC_FORMAT)
# handler.setFormatter(formatter)
# root = logging.getLogger(__name__)
# root.setLevel(os.environ.get("LOGLEVEL", "INFO"))
# root.addHandler(handler)

# logging.basicConfig(filename=os.path.join(os.getcwd(),
#                                           os.path.basename(__file__).rstrip('.py'),
#                                           '.log'), filemode='w')


def main(logoutto: str):
    fmtstr = " Name: %(user_name)s : %(asctime)s: (%(filename)s): %(levelname)s: %(funcName)s Line: %(lineno)d - %(message)s"
    datestr = "%m/%d/%Y %I:%M:%S %p "

    # basic logging config
    logging.basicConfig(
        filename=os.path.join(os.getcwd(),
                              os.path.basename(logoutto).rstrip('.py') + '.log'),
        level=logging.DEBUG,
        filemode="w",
        format=fmtstr,
        datefmt=datestr,
    )


def stdout(logoutto: str):
    sys.stdout = open(os.path.join(os.getcwd(),
                              os.path.basename(logoutto).rstrip('.py') + '.out'),
                      'w')
