import sys
import os.path
from math import log

def approx_equal(expected, observed, tolerance=0.05):
    v = abs(log(observed / float(expected)))
    return v <= tolerance

def add_root_path():
    'force phenoseq source root to be included in sys.path'
    p = os.path.dirname(os.path.abspath(__file__)) # phenoseq test dir
    sys.path.append(os.path.dirname(p)) # add phenoseq source root dir

