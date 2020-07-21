import sys
import pickle
import inspect


class MyPrintClass(object):

    def __init__(self, debug_mode):
        self.debug = debug_mode

    def myprint(self, *objects, sep=' ', end='\n', file=sys.stdout, flush=False):
        if self.debug:
            print(*objects, sep=sep, end=end, file=file, flush=flush)

