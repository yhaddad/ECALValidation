__author__  = 'Yacine Haddad'
__email__   = 'yhaddad@cern.ch'
__version__ = '0.1.0'


import matplotlib.pyplot as plt
import pandas     as pd
import numpy      as np
import ParseTable as pt
import ROOT       as r
import matplotlib.dates as dates

from pprint       import pprint
from root_numpy   import root2array, tree2array, hist2array
from rootpy.io    import root_open
from collections  import OrderedDict
from termcolor    import colored
from jsmin        import jsmin

import peakutils
import os, sys, json,logging

from scipy.interpolate   import UnivariateSpline, InterpolatedUnivariateSpline, LSQUnivariateSpline
from scipy               import signal, stats
from scipy.interpolate   import splev, splrep
from   matplotlib.ticker import NullFormatter
from   matplotlib.ticker import MultipleLocator, FormatStrFormatter
