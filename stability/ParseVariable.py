from __init__ import *

logging.basicConfig(format=colored('%(levelname)s:', attrs=['bold'])
                    + colored('%(name)s:', 'blue') + ' %(message)s')
logger = logging.getLogger('MoECAL')
logger.setLevel(level=logging.INFO)

class variable():
    def __init__(self, name="", options={}):
        self.__template__ = {
            "name"     : "",
            "range"    : "",
            "title"    : "",
            "bins"     : 100,
            "draw"     : True,
            "adjust"   : False
        }
        self.__dict__ = self.__template__
        self.__dict__.update(options)
        self.name = name
    def __str__(self):
        if self.draw:
            return colored(" -- variable :: %18s -- %12s" % (self.name, self.draw ), "green" )
        else:
            return colored(" -- variable :: %18s -- %12s" % (self.name, self.draw ), "red" )

class monitor():
    def __init__(self,
                 config="config/basic_ecal_configuration.json",
                 run_ranges="data-stability/80X_dataRun2_Prompt_v11_Nov2016_interval_100000.dat"):
        assert os.path.exists(run_ranges)==True, colored("[ERROR] ran range file doesn't exist ... ","red")
        self.variables      = {}
        self.regions        = {}
        self.ecal_selection = {}
        self.ecal_regions   = {}
        self.columns        = []
        self.data           = {}
        self.run_ranges     = pt.read_run_range(path=os.path.dirname (run_ranges),
                                                file=os.path.basename(run_ranges))
        self._read_ecal_configuration_(infile=config)
        self.columns = ['eventNumber', 'runNumber', 'lumiBlock'] +  self.variables.keys()
        self.monitor_features = []
        print colored("------- run-ranges :", "cyan" , "on_white")
        print self.run_ranges.run_number

    def _read_ecal_configuration_(self, infile):
        assert os.path.exists(infile)==True, colored("[ERROR] configuration file doesn't exist ... ","red")
        _config_ = None
        with open(infile) as f:
            _config_ = json.loads(f.read())
        for key in _config_:
            if 'variables' in key.lower():
                logger.info(colored("------- variables :", "yellow"))
                for var in _config_[key]:
                    v = variable(var, _config_[key][var])
                    self.variables[v.name] = v
                    logger.info(v)
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
        _data_ = { region : None for region in self.ecal_regions}
        logger.info(colored("------- ntuples :", "yellow"))
        for index, root in _config_.iterrows():
            chain = r.TChain('merged_' + str(index) )
            chain.Add(root.file+'/'+root.tree)
            _df_ = None
            for region,cut in self.ecal_regions.items():
                _cut_ = "&&".join( [cut, self.ecal_selection[selection]] )
                _df_  = tree2array( chain, selection = _cut_ , branches = self.columns )
                logger.info(colored(" -- ntuple :: %10s -- %10s -- %12i" % ( root.id, region ,_df_.shape[0] ), "green" ))
                if _data_.get(region, None) == None:
                    _data_[region] = _df_
                else:
                    _data_[region] = np.concatenate((_data_[region],  _df_), axis=0)
        self.data = self._flatten_data_(_data_)
        return self.data


    def _flatten_data_(self, indata = None):
        assert indata is not None, logger.error(colored("[ERROR] the input data seems to be None .... ","yellow","on_red"))
        _data_  = {}
        for region, dd in indata.items():
            df = pd.DataFrame()
            for name in dd.dtype.fields:
                if name not in self.columns : continue
                if len(dd[name].shape) > 1 :
                    for i in range(dd[name].shape[1]):
                        df[name + '_' + str(i)] = dd[name][:,i]
                        pprint(self.variables[name].__dict__)
                        self.variables[name + '_' + str(i) ] = variable(name, self.variables[name].__dict__.copy())
                        self.variables[name + '_' + str(i) ].name = name + '_' + str(i)
                        self.variables[name + '_' + str(i) ].draw = True
                        self.variables[name].draw = False
                else:
                    df[name] = dd[name]
            _data_[region] =  df
        logger.info(colored("------- variables status:", "yellow" ))
        for key, var in self.variables.items() :
            logger.info( var )
        return _data_

    def fit(self, draw_fit=True, outdir = './'):
        if not os.path.exists(outdir): os.makedirs(outdir)
        fig = plt.figure(figsize=(5,7))
        for region in self.ecal_regions:
            print colored("------- drawing : %s" % region , "yellow" )
            for v, var in self.variables.items():
                logger.info('\t' + var.__str__() )
                if var.draw == False : continue
                for index, row in self.run_ranges.iterrows():
                    _data_ = self.data[region][np.logical_and(self.data[region].runNumber >= row.run_min,
                                                              self.data[region].runNumber <= row.run_max,
                                                             )]
                    if not os.path.exists(os.path.join(outdir,v,'fits')): os.makedirs(os.path.join(outdir,v,'fits'))
                    print colored("\t --  : %15s %3i %10i %10s" % (v, index, _data_.shape[0], row.run_number ) , "green" )

                    if _data_.shape[0] != 0 :
                        y,xbin = np.histogram(_data_[v], bins=var.bins, range=var.range)
                        if y.sum() == 0 :
                            print colored('[warning] :: the variable [%s] is empty !! ' % v, 'red')
                            var.draw = False
                            continue
                        x    = [(xbin[i]+xbin[i+1])/2.0   for i in range(0,len(xbin)-1)]
                        yerr = np.array([np.sqrt(y[i]) for i in range(0,len(y) )])
                        yerr[yerr == 0] = np.sqrt(-np.log(0.68))
                        peak,pmin,pmax,_ = self.find_FWHM(x,y,yerr, xrange=var.range, draw=True,
                                                          label='spline-%s-%i-%i-%s'    % (v, row.run_min, row.run_max, region ),
                                                          title='%s run-range = (%i-%i)'% (region, row.run_min, row.run_max))

                        plt.savefig(os.path.join(os.path.join(outdir,v,'fits'),
                                     'variable-%s-%i-%i-%s.png' % (v, row.run_min, row.run_max, region) ))
                        plt.clf()
                        self.run_ranges.loc[index, '%s_%s_mean' % (v,region)] = _data_[v].mean()
                        self.run_ranges.loc[index, '%s_%s_std'  % (v,region)] = _data_[v].std()
                        self.run_ranges.loc[index, '%s_%s_krts' % (v,region)] = stats.kurtosis(_data_[v])
                        self.run_ranges.loc[index, '%s_%s_skew' % (v,region)] = stats.skew(_data_[v])
                        self.run_ranges.loc[index, '%s_%s_peak' % (v,region)] = peak
                        self.run_ranges.loc[index, '%s_%s_pmin' % (v,region)] = pmin
                        self.run_ranges.loc[index, '%s_%s_pmax' % (v,region)] = pmax
                    else:
                        self.run_ranges.loc[index, '%s_%s_mean' % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_std'  % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_krts' % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_skew' % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_peak' % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_pmin' % (v,region)] = -999.0
                        self.run_ranges.loc[index, '%s_%s_pmax' % (v,region)] = -999.0
        print os.path.join(outdir, 'fit_data.tex')
        fout_ = open(os.path.join(outdir, 'fit_data.tex'), 'w')
        fout_.write(self.run_ranges.to_latex())
        fout_.close()


    def monitor_peak(self, var=None, region=None, outdir='./plots'):
        datavalues = self.run_ranges['%s_%s_peak' % (var.name, region)]
        xdata      = self.run_ranges.run_number

        left, width    = 0.1, 1.0
        bottom, height = 0.1, 0.5
        hist_sep       = 0.0
        rect_hist = [left+width+hist_sep, bottom, 0.25*height, height]
        rect_plot = [left, bottom, width, height]

        nullfmt = NullFormatter()
        fig = plt.figure(figsize=(12,6))

        ax_plot = plt.axes(rect_plot)
        ax_hist = plt.axes(rect_hist)
        for k, spine in ax_plot.spines.items():
            spine.set_zorder(10)
        for k, spine in ax_hist.spines.items():
            spine.set_zorder(10)

        ax_plot.set_title ( 'monitoring : %s %s' % (region, var.name) )
        ax_plot.set_ylabel( var.title )
        ax_hist.yaxis.set_major_formatter(nullfmt)
        y,_,_ = ax_hist.hist(datavalues, bins=var.bins,
                             orientation = 'horizontal' ,
                             histtype    = 'stepfilled' , alpha=0.6,zorder=10,)
        xbins = range(1,1+len(datavalues),1)
        xerrs = 0.5 * np.ones(self.run_ranges.shape[0])
        ax_plot.errorbar(xbins,datavalues,color='blue',xerr = xerrs,label='peak',
                         capthick=0,marker='.',ms=6,ls='None',zorder=10,)
        err_ = (self.run_ranges[ '%s_%s_pmax' % (var.name,region) ] - self.run_ranges[ '%s_%s_pmin' % (var.name,region) ])
        pos_ = (self.run_ranges[ '%s_%s_pmax' % (var.name,region) ] + self.run_ranges[ '%s_%s_pmin' % (var.name,region) ])/2.0

        ax_plot.bar( xbins ,   err_,
                bottom = pos_ - err_ / 2.0,
                width  = 2*xerrs ,
                color='blue'  ,alpha=0.3, zorder=9,label='FWHM',
                align='center',edgecolor='None',lw=0)

        majorLocator = MultipleLocator(2)
        minorLocator = MultipleLocator(1)
        ax_plot.xaxis.set_major_locator(majorLocator)
        ax_plot.xaxis.set_minor_locator(minorLocator)
        xlabels = ax_plot.get_xticks().tolist()
        for i in range(2,len(xdata)+2,2):
            xlabels[i/2+1] = xdata.tolist()[i-1]
        for i in range(len(xlabels)):
            if type(xlabels[i]) == type(1.0): xlabels[i] = ''
        ax_plot.set_xticklabels(xlabels, family="monospace")
        xlabels = ax_plot.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=8)

        ax_plot.xaxis.grid(True, which="minor")
        ax_plot.yaxis.grid()
        ax_hist.xaxis.grid(True, which="minor")
        ax_hist.yaxis.grid()
        ax_hist.xaxis.set_ticks([])

        if var.adjust==False :
            ax_plot.set_ylim(var.range)
            ax_hist.set_ylim(var.range)

        legend = ax_plot.legend(loc='best')

        ax_hist.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
        ax_plot.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
        ax_plot.grid(which='minor', color='0.85', linestyle='--',dashes=(5,1),zorder=0)
        outdir = os.path.join(outdir,var.name)
        if not os.path.exists(outdir): os.makedirs(outdir)

        plt.savefig(os.path.join(outdir, '%s_%s_peak.png' % (var.name, region)), format='png',orientation='landscape',
                    dpi=200,papertype='a4',pad_inches=0.1,
                    bbox_inches='tight')
        plt.savefig(os.path.join(outdir, '%s_%s_peak.pdf' % (var.name, region)), format='png',orientation='landscape',
                    dpi=200,papertype='a4',pad_inches=0.1,
                    bbox_inches='tight')
    def monitor_mean(self, var=None, region=None, outdir='./plots'):
            datavalues = self.run_ranges['%s_%s_mean' % (var.name, region)]
            xdata      = self.run_ranges.run_number

            left, width    = 0.1, 1.0
            bottom, height = 0.1, 0.5
            hist_sep       = 0.0
            rect_hist = [left+width+hist_sep, bottom, 0.25*height, height]
            rect_plot = [left, bottom, width, height]

            nullfmt = NullFormatter()

            fig = plt.figure(figsize=(12,6))

            ax_plot = plt.axes(rect_plot)
            ax_hist = plt.axes(rect_hist)
            for k, spine in ax_plot.spines.items():
                spine.set_zorder(10)
            for k, spine in ax_hist.spines.items():
                spine.set_zorder(10)
            ax_plot.set_title ( r'monitoring mean: %s %s' % (region, var.name) )
            ax_plot.set_ylabel( r'$\langle %s \rangle$' % var.title.replace('$',r'') )
            ax_hist.yaxis.set_major_formatter(nullfmt)
            y,_,_ = ax_hist.hist(datavalues, bins=var.bins,
                                 orientation ='horizontal',
                                 histtype    ='stepfilled', alpha=0.6,zorder=10,)
            xbins = range(1,1+len(datavalues),1)
            xerrs = 0.5 * np.ones(self.run_ranges.shape[0])
            ax_plot.errorbar(xbins,datavalues,color='blue',xerr = xerrs,
                             label='peak',yerr=np.std(datavalues)/np.sqrt(datavalues),
                             capthick=0,marker='.',ms=6,ls='None',zorder=10,)
            majorLocator = MultipleLocator(2)
            minorLocator = MultipleLocator(1)
            ax_plot.xaxis.set_major_locator(majorLocator)
            ax_plot.xaxis.set_minor_locator(minorLocator)
            xlabels = ax_plot.get_xticks().tolist()
            for i in range(2,len(xdata)+2,2):
                xlabels[i/2+1] = xdata.tolist()[i-1]
            for i in range(len(xlabels)):
                if type(xlabels[i]) == type(1.0): xlabels[i] = ''

            ax_plot.set_xticklabels(xlabels, family="monospace")
            xlabels = ax_plot.get_xticklabels()
            plt.setp(xlabels, rotation=90, fontsize=8)

            ax_plot.xaxis.grid(True, which="minor")
            ax_plot.yaxis.grid()
            ax_hist.xaxis.grid(True, which="minor")
            ax_hist.yaxis.grid()
            ax_hist.xaxis.set_ticks([])
            if var.adjust==False :
                ax_plot.set_ylim(var.range)
                ax_hist.set_ylim(var.range)

            legend = ax_plot.legend(loc='best')

            ax_hist.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
            ax_plot.grid(which='major', color='0.7' , linestyle='--',dashes=(5,1),zorder=0)
            ax_plot.grid(which='minor', color='0.85', linestyle='--',dashes=(5,1),zorder=0)

            outdir = os.path.join(outdir,var.name)
            if not os.path.exists(outdir): os.makedirs(outdir)

            plt.savefig(os.path.join(outdir, '%s_%s_mean.png' % (var.name, region)), format='png',orientation='landscape',
                        dpi=200,papertype='a4',pad_inches=0.1,
                        bbox_inches='tight')
            plt.savefig(os.path.join(outdir, '%s_%s_mean.pdf' % (var.name, region)), format='png',orientation='landscape',
                        dpi=200,papertype='a4',pad_inches=0.1,
                        bbox_inches='tight')

    def find_FWHM(self, x,y, yerr, xrange, draw=True, label= '', title=''):
        t = np.linspace(xrange[0],xrange[1],1000)
        w = 1.0/yerr
        w[w == np.inf] = 0
        spl   = UnivariateSpline(x, y,w=w)
        maximum = spl(x)[np.argmax(spl(x))]
        id0     = peakutils.indexes(spl(t), thres=1/max(spl(t)), min_dist=50)
        # find the width
        spline_  = UnivariateSpline(x,y-maximum/2,w=w)
        roots_   = spline_.roots()
        peak_    = t[np.argmax(spline_(t))]
        peaks_   = spl.roots()
        p0       = [t[p] for p in id0]
        if draw:
            plt.subplots_adjust(hspace=0)
            ax1 = plt.subplot2grid((4,1), (0,0),rowspan=2)
            ax1.set_ylabel(r"events")
            # logger.debug("\t UnivariateSpline : chi2/ndf == %1.3f " % spl.get_residual()/len(y))
            # logger.debug("\t UnivariateSpline : Ncoef    == %i " % len(spl.get_coeffs()))
            # logger.debug("\t UnivariateSpline : Nknots   == %i " % len(spl.get_knots()))
            # logger.debug("\t UnivariateSpline : smooth   == %i " % len(w))
            ax1.fill_between(t,spl(t)-np.sqrt(spl(t)),spl(t)+np.sqrt(spl(t)),
                             alpha=0.5, edgecolor='#f2ae72', facecolor='#f2ae72')

            ax1.errorbar(x,y,yerr=yerr, capthick=0,marker='.', ms=4,ls='None',c='Blue')
            ax1.plot(t,spl(t) ,'r-' ,alpha=0.9)
            for p in roots_: ax1.axvline(x=p)
            ax1.axvline(x=peak_,color='red')
            plt.ylim([0,max(y)*1.2])
            ax2 = plt.subplot2grid((4,1), (2,0),sharex=ax1)
            # ----------------------------
            residuals = (spl(x)-y)/spl(x)
            res_min   = -(np.sqrt(spl(t))/spl(t))
            res_max   = +(np.sqrt(spl(t))/spl(t))
            ax2.fill_between(t,res_min,res_max, alpha=0.5, edgecolor='#f2ae72', facecolor='#f2ae72')

            ax2.errorbar(x, residuals, yerr=yerr/spl(x),capthick=0,marker='.',
                         ms=4,ls='None',c='Blue')
            plt.ylim([-1.5,1.5])

            ax2.axhline(y=0,color='Red')
            ax3 = plt.subplot2grid((4,1), (3,0),sharex=ax2)
            cum_y   = np.cumsum(y)       / float(np.sum(y))
            cum_spl = np.cumsum(spl(t))/ np.sum(spl(t))
            ax3.plot(x,cum_y   ,'b-' ,alpha=0.9)
            ax3.plot(t,cum_spl ,'r-' ,alpha=0.9)
            Dn = np.max(np.abs(cum_y - (np.cumsum(spl(x))/ np.sum(spl(x)))))
            ax3.text(0.3,0.9,'KS Test : %1.4f' % Dn, horizontalalignment='right',
                     verticalalignment='top',
                     transform=ax3.transAxes,fontsize=7)
            plt.ylim([0,1])
            ax1.xaxis.set_visible(False)
            ax2.xaxis.set_visible(False)

            plt.xlim(xrange)
            plt.show()
        return peak_ , min(roots_), max(roots_), spl.get_residual()/len(y)
