import matplotlib.pyplot as plt
import pandas     as pd
import numpy      as np
import ParseTable as pt
import ROOT       as r
from pprint import pprint
from root_numpy import root2array, tree2array
from rootpy.io  import root_open
from scipy.interpolate import UnivariateSpline
from collections  import OrderedDict
from termcolor    import colored
from jsmin        import jsmin
import os, sys, json

class variable():
    def __init__(self, name="", options={}):
        self.__template__ = {
            "name"     : "",
            "range"    : "",
            "title"    : "",
            "bins"     : 100
        }
        self.__dict__ = self.__template__
        self.__dict__.update(options)
        self.name = name
    def __str__(self):
        return colored(" -- variable :: %18s -- %12s" % (self.name, "monitor" ), "green" )

class monitor():
    def __init__(self, config="config/basic_ecal_configuration.json"):
        self.variables      = {}
        self.regions        = {}
        self.ecal_selection = {}
        self.ecal_regions   = {}
        self._read_ecal_configuration_(infile=config)
        self.columns = []
        self.columns = ['eventNumber', 'runNumber', 'lumiBlock'] +  self.variables.keys()

    def _read_ecal_configuration_(self, infile):
        assert os.path.exists(infile)==True, "[ERROR] path error, or file doesn't exist ... "
        _config_ = None
        with open(infile) as f:
            _config_ = json.loads(f.read())#jsmin(f.read()))
        for key in _config_:
            if 'variables' in key.lower():
                print colored("------- variables :", "green" )
                for var in _config_[key]:
                    v = variable(var, _config_[key][var])
                    self.variables[v.name] = v
                    print v
            if "selection" in key.lower():
                self.ecal_selection = _config_[key]
            if "regions"    in key.lower():
                self.ecal_regions   = _config_[key]

    def read_ntuple(self, path="./", cfg="", selection = 'loose'):
        """
        read the config files the same as for the Zfitter
        to the a chain of the files
        """
        _config_ =  pd.read_csv(path + "/" + cfg , sep = "\t", names = ['id', 'tree', 'file'], comment ="#")
        chain = r.TChain('merged')
        print colored("------- ntuples :", "green" )
        for index, root in _config_.iterrows():
            print colored(" -- ntuple :: %10s -- %12s" % ( root.id, root.file ), "green" )
            chain.Add(root.file+'/'+root.tree)
        _data_ = {}
        for region,cut in self.ecal_regions.items():
            _cut_ = "&&".join( [cut, self.ecal_selection[selection]] )
            _data_[region] = tree2array( chain, selection = _cut_ , branches = self.columns )
        dd = self._flatten_data_(_data_)
        return dd

    def _flatten_data_(self, indata = None):
        assert indata is not None, colored("[ERROR] the input data seems to be None .... ","red")
        _data_ = {}
        for region, dd in indata.items():
            df = pd.DataFrame()
            for name in dd.dtype.fields:
                if name not in self.columns : continue
                if len(dd[name].shape) > 1 :
                    for i in range(dd[name].shape[1]):
                        df[name + '_' + str(i)] = dd[name][:,i]
                else:
                    df[name] = dd[name]
            _data_[region] =  df
        return _data_
    def save_fit_plots(self, outdir = './'):
        pass

    


# def run_test():
#     run_ranges = pt.read_run_range(path='./ntuples/',file='ICHEP_interval_100000.dat')
#     print run_ranges.head()
#     raw_data = read_ntuple('./', 'config.dat', 'loose25nsRun2')
#     print pprint(raw_data)
#     data = flatten_data(indata = raw_data)
#     for reg in regions :
#         print "region :: ", reg
#         pprint(data[reg].head())
