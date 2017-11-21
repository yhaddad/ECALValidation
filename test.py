from stability import ParseVariable as pv

mm = pv.monitor(
    config='config/basic_ecal_configuration.json',
    run_ranges='./monitoring_HighEta_interval.dat',
)

mm.read_ntuple(
    path="./",
    cfg="monitoring_2017_Z_highEta.dat",
    selection='simple'
)

# print dir(mm)
# print type(mm.data)

data = mm.run_ranges
# for i, c in mm.data.items():
#     print i
#     c.to_hdf('raw_data_92X_%s.h5' % i,  'monitoring', mode='w', format='table')

# for i, c in mm.simu.items():
#     print i
#     c.to_hdf('raw_simu_92X_%s.h5' % i,  'monitoring', mode='w', format='table')

mm.fit(outdir='92X_dataRun2_Prompt_v11')
data.to_hdf('data_92X_dataRun2_Prompt_v11.h5',
            'monitoring', mode='w', format='table')

for v, var in mm.variables.items():
    for region in mm.ecal_regions:
       mm.monitor_peak(var, region, outdir='./tmp')
       mm.monitor_mean(var, region, outdir='./tmp') 

"""
def bwias_peak_monitor(xbin, hist, x, spl, label='', args = None):
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
    ax1.text(0.95,0.90,'peak = {0:1.2f} GeV'.format((peaks[0])),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax1.transAxes,fontsize=7)
    ax1.text(0.95,0.85,'FWHM = {0:1.2f} GeV'.format((max(roots)-min(roots))),
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
    plt.toring_2017_Z_golden_interval_10000.datshow()
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
"""
