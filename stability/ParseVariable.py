import matplotlib.pyplot as plt
import pandas     as pd
import numpy      as np
import ParseTable as pt
import ROOT       as r
from root_numpy import root2array, tree2array
from rootpy.io  import root_open

range_header     = ["run","nevent","date"]
ecal_selections  = {
    "loose25nsRun2" : "(((eleID[0] & 131072)==131072)&&((eleID[1] & 131072)==131072))&&((energySCEle_must_regrCorr_ele[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_must_regrCorr_ele[1]/cosh(etaSCEle[1]) >= 25))",
    "loose" : "(((eleID[0] & 128)==128)&&((eleID[1] & 128)==128))&&((energySCEle_must_regrCorr_ele[0]/cosh(etaSCEle[0]) >= 25)&&(energySCEle_must_regrCorr_ele[1]/cosh(etaSCEle[1]) >= 25))"
}

regions = {
    "EB"   : "abs(etaEle_ele[0]) < 1.4442 && abs(etaEle_ele[1]) < 1.4442",
    "EE"   : "abs(etaEle_ele[0]) > 1.4442 && abs(etaEle_ele[1]) > 1.4442",
}

variables = {
    "nPV"           : {"range":[0,30 ] },
    "seedXSCEle[0]" : {"range":[0,500] },
    "seedXSCEle[1]" : {"range":[0,500] }
}


data = None

def read_ntuple(path="./", cfg="", selection = 'loose'):
    # read the config files the same as for the Zfitter
    # to the a chain of the files
    config =  pd.read_csv(path + "/" + cfg , sep = " ", names = ['id', 'tree', 'file'], comment ="#")
    print config
    chain = r.TChain('merged')
    for index, root in config.iterrows():
        print index, "  merging the file : ", root.file
        chain.Add(root.file+'/'+root.tree)
    # transform this chain to an array ment to be used later by matplotlib
    data = tree2array( chain,
                       selection =  ecal_selections[selection],
                      )
    return data

def monitor_variable(path, file, variable):
    run_ranges = pt.read_run_range(path=path,file=file)

    run_ranges[variable + '_mean'] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_std' ] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_min' ] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_max' ] = np.zeros(run_ranges.shape[0])

    for index, row in run_ranges.iterrows():
        run_ranges['_mean'][index] = data[np.logical_and(data['runNumber']>=row.run_min, data['runNumber']<=row.run_max)]['nPV'].mean()
        run_ranges['_std' ][index] = data[np.logical_and(data['runNumber']>=row.run_min, data['runNumber']<=row.run_max)]['nPV'].std()
        run_ranges['_std' ][index] = data[np.logical_and(data['runNumber']>=row.run_min, data['runNumber']<=row.run_max)]['nPV'].size()
        print run_ranges['run_number'][index], ' -- ', run_ranges['nPV_mean'][index]

    run_ranges['_sem'  ] = run_ranges['_mean']/np.sqrt(run_ranges['Nevents'])
    run_ranges['_stde' ] = run_ranges['_std' ]/np.sqrt(run_ranges['Nevents'])

    run_ranges['bin']       = range(0,run_ranges.shape[0])
    run_ranges['bin_error'] = 0.5 * np.ones(run_ranges.shape[0])

    for k, spine in ax.spines.items():
        spine.set_zorder(10)
    ax.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
    ax.grid(which='minor', color='0.85', linestyle='--',dashes=(5,1),zorder=0)

    ax.errorbar(run_ranges['bin'],run_ranges['nPV_mean'],
             yerr=run_ranges['_sem'  ],
             xerr=run_ranges['_error'],
             capthick=0,marker='o',ms=4,ls='None', zorder=10)
    plt.ylim([0,60])
    plt.savefig('test.pdf')


def test():
    data = read_ntuple(path="./", cfg="config.dat", selection = 'loose')
    monitor_variable(path='./ntuples/',file='ICHEP_interval_100000.dat', variable='nPV')
