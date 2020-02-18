"""
Automatically adds the project's files to the PYTHONPATH
"""

import os, sys

this_dir = os.path.dirname(__file__)
home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '1000Genomes')

if proj_dir not in sys.path:
  sys.path.insert(0, proj_dir)