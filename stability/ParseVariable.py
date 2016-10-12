import matplotlib.pyplot as plt
import pandas     as pd
import numpy      as np
import ParseTable as pt
import ROOT       as r
from root_numpy import root2array, tree2array
from rootpy.io  import root_open
from scipy.interpolate import UnivariateSpline


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


df = None

def read_ntuple(path="./", cfg="", selection = 'loose'):
    # read the config files the same as for the Zfitter
    # to the a chain of the files
    config =  pd.read_csv(path + "/" + cfg , sep = " ", names = ['id', 'tree', 'file'], comment ="#")
    print config
    chain = r.TChain('merged')
    for index, root in config.iterrows():
        print root.id , "\t: ", root.file
        chain.Add(root.file+'/'+root.tree)
    # transform this chain to an array ment to be used later by matplotlib

    data = tree2array( chain,
                       selection =  ecal_selections[selection],
                      )
    return data
def find_FWHM(bins, hist):
    maximum = hist[np.argmax(hist)]
    peak    = bins[np.argmax(hist)]
    spline  = UnivariateSpline(bins[:-1], hist-maximum/2 , s=0)
    return peak, spline.roots()

def bias_peak_monitor(xbin, hist, x, spl, label='', args = None):
    if args is None:
        args = {}
    fig = plt.figure(figsize=(4,4.5))
    plt.subplots_adjust(hspace=0)
    # ----------------------------------
    maximum = args.get('maximum', max(hist))
    peaks   = args.get('peaks'  ,[np.mean(hist)])
    roots   = args.get('roots'  ,[0,10])
    # ----------------------------------
    ax1 = plt.subplot2grid((3,1), (0,0),rowspan=2)
    ax1.text(0.95,0.95,'Range : ' + label,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax1.transAxes,fontsize=7)
    ax1.text(0.95,0.90,'peak = %1.2f GeV'% (peaks[0]),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax1.transAxes,fontsize=7)
    ax1.text(0.95,0.85,'FWHM = %1.2f GeV'% (max(roots)-min(roots)),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax1.transAxes,fontsize=7)
    ax1.plot(xbin[:-1], hist,'.')
    ax1.plot(x        , spl(x)  )
    for p in peaks: ax1.plot((p,p), (0,maximum), 'r-')
    ax1.axvspan(min(roots),max(roots), facecolor='b', alpha=0.2)
    plt.ylim([0,maximum*1.2])

    ax1.tick_params(
        axis        = 'x'   , # changes apply to the x-axis
        which       = 'both', # both major and minor ticks are affected
        bottom      = 'off' , # ticks along the bottom edge are off
        top         = 'off' , # ticks along the top edge are off
        labelbottom = 'off' ) # labels along the bottom edge are off
    ax2 = plt.subplot2grid((3,1), (2,0),sharex=ax1)
    ax2.axvspan(min(roots),max(roots), facecolor='b', alpha=0.2)
    # ----------------------------
    residuals = (spl(xbin[:-1])-hist)/spl(xbin[:-1])
    ax2.plot(xbin[:-1], residuals , 'r.' )
    plt.ylim([-1.5,1.5])
    ax2.axhline(y=0)
    ax2.axhspan(-0.5,0.5, facecolor='b', alpha=0.2)
    plt.show()
    plt.savefig('log-plots/univariate-spline-' + label + '.png')
    plt.savefig('log-plots/univariate-spline-' + label + '.pdf')


def monitor_variable(df, path, file, variable):
    run_ranges = pt.read_run_range(path=path,file=file)

    run_ranges[variable + '_mean'] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_std' ] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_peak'] = np.zeros(run_ranges.shape[0])
    run_ranges[variable + '_nevt'] = np.zeros(run_ranges.shape[0])
    for index, row in run_ranges.iterrows():
        _data_ = df[np.logical_and(df['runNumber']>=row.run_min, df['runNumber']<=row.run_max)]['nPV']
        run_ranges[variable + '_mean'][index] = _data_.mean()
        run_ranges[variable + '_std' ][index] = _data_.std()
        run_ranges[variable + '_nevt'][index] = _data_.size
        # find a maximum
        bins,hist = np.histogram(_data_,n)

        print run_ranges['run_number'][index], ' -- ', run_ranges['nPV_mean'][index]

    run_ranges[variable + '_sem'  ] = run_ranges[variable + '_mean']/np.sqrt()
    run_ranges[variable + '_stde' ] = run_ranges[variable + '_std' ]/np.sqrt()

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
    data = read_ntuple(path="./", cfg="config.dat", selection = 'loose25nsRun2')
    monitor_variable(df=data, path='./ntuples/',file='ICHEP_interval_100000.dat', variable='nPV')
