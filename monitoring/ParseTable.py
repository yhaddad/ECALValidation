# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.ticker import NullFormatter
from   matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.dates as dates
import pandas as pd
import glob, os, re
import datetime
import re
import string
from math import sqrt

region_labels = {
    "EB"                 : "EB",
    "EB-gold"            : "EB $R_{9} > 0.94$",
    "EB-bad"             : "EB $R_{9} < 0.94$",
    "EB-absEta_0_1"      : "EB $\mid\eta\mid < 1$",
    "EB-absEta_1_1.4442" : "EB $1 < \mid\eta\mid < 1.4442$",
    "EE"      : "EE",
    "EE-gold" : "EE $R_{9} > 0.94$",
    "EE-bad" : "EE $R_{9} < 0.94$",
    "EE-absEta_1.566_2" : "EE $1.566 < \mid\eta\mid < 2$",
    "EE-absEta_2_2.5" : "EE $2 < \mid\eta\mid < 2.5$",
    "inclusive" : "All $\mid\eta\mid$ and $R_{9}$",
    "gold" : "All $\mid\eta\mid$  $R_{9} > 0.94$",
    "bad" : "All $\mid\eta\mid$ $R_{9} < 0.94$"
}

var_labels = {
    "DeltaM_data"        : "$\Delta m$",
    "DeltaP"             : "$\Delta p$",
    "width_data"         : "$\sigma_{CB}$",
    "rescaledWidth_data" : "$\sigma_{CB}$ (Rescaled)",
    "additionalSmearing" : "Additional Smearing",
    "chi2data"           : "$\chi^{2}$",
    "events_lumi"        : "events/lumi",
    "sigmaeff_data"      : "$\sigma_{eff}$",
    "sigmaeff_data_thirty"      : "$\sigma_{eff30}$",
    "sigmaeff_data_fifty"      : "$\sigma_{eff50}$"
}

format_fig_output = ['pdf','png']

def read_regions_from_table(path = "",tableFile= "",xVar=""):

    regions = []
    with open(path+'/'+tableFile) as f:
        for line in f.read().split('\n'):
            if line == '': continue
            if line[0] == '#': continue
            if 'category' in line: continue

            region = line.split('&')[0]
            if xVar == '':
                if region.split('-runNumber')[0] not in regions:
                    regions.append(region.split('-runNumber')[0])
            else:
                part = re.findall('[%s]+' % string.ascii_letters,region.replace(xVar,''))
                if len(part) == 0:
                    part = 'inclusive'
                    if part not in regions:
                        regions.append(part)
                else:
                    if part[0] not in regions:
                        regions.append(part[0])

    return regions


def read_regions_from_regionsfile(path = "",regionsFile=""):

    assert regionsFile != "", "The regions file must be specified"

    regions = []
    with open(path+regionsFile) as f:
        for line in f.read().split('\n'):
            if line != '' and line[0] != '#':
                regions.append(line.split()[0])

    return regions


def parse_table_over_regions(path = "", tableFile = "",category="",xVar=""):

    assert xVar != "", "The x variable name must be specified!"

    if category == "inclusive":
        cats = ['EB','EE','gold','bad'] #if inclusive specify the general eta and R9 regions to check against
    else:
        cats = [category]

    lines = []
    varmins = []
    varmaxes = []

    text = open(path+tableFile).read()
    text = re.sub('\\pm','~',text)
    text = re.sub(r'[\\£@#${}^]','',text)

    assert xVar in text, "The specified variable label must be in the stability file!"

    #Parse top line for variables
    variables = text.split("\n")[0].split(' & ')

    for line in text.split('\n'):

        if line == '': continue
        if line[0] == '#': continue
        if 'category' in line: continue

        if category == 'inclusive':
            if not any(c in line.split('&')[0] for c in cats):
                lines.append(line)
                varmins.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))
                varmaxes.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]))
        else:
            if category in line.split('&')[0]:
                lines.append(line)
                varmins.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[0]))
                varmaxes.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",line.split('&')[0])[1]))

    #If there's nothing for this category, return an empty dataframe
    if len(lines) == 0: return pd.DataFrame()

    data = pd.DataFrame()
    data[xVar + '_min'] = varmins
    data[xVar + '_max'] = varmaxes

    variables = [xVar+'min',xVar+'max'] + variables[2:]

    for i in range(2,len(variables),1):
        values = []
        errors = []
        for line in lines:
            value = line.split('&')[i].replace(' ','')

            if '--' in value: value = '0'

            if len(value.split('~')) == 2:
                values.append(float(value.split('~')[0]))
                errors.append(float(value.split('~')[1]))
            elif len(value.split('~')) == 1:
                values.append(float(value))
                errors.append(0.0)

            if 'mc' in variables[i]:
                variables[i] = variables[i].replace("mc","MC")
            if '/' in variables[i]:
                variables[i] = variables[i].replace("/","_")

        data[variables[i]] = values
        data[variables[i]+'_err'] = errors

    return data


def read_run_range(path = "", file = ""):

    raw_    = {
        'run_number': [],
        'Nevents'   : [],
        'UnixTime'  : []
    }

    with open(path+file) as f:
        for line in f.read().split('\n'):
            if len(line)==0: continue
            values = line.split('\t')
            raw_['run_number'].append(values[0])
            raw_['Nevents'   ].append(values[1])
            raw_['UnixTime'  ].append(values[2])
    data_ = pd.DataFrame(raw_)

    #Transform the entries
    data_.Nevents = [float(x) for x in data_.Nevents]
    data_[ "UnixTime_min" ] = [ int(x.split('-')[0]) for x in data_.UnixTime ]
    data_[ "UnixTime_max" ] = [ int(x.split('-')[1]) for x in data_.UnixTime ]
    data_[ "run_min" ] = [ int(x.split('-')[0]) for x in data_.run_number ]
    data_[ "run_max" ] = [ int(x.split('-')[1]) for x in data_.run_number ]
    data_[ "date_min"] = [ datetime.datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in data_.UnixTime_min ]
    data_[ "date_max"] = [ datetime.datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%M:%S') for x in data_.UnixTime_max ]
    data_['time'] = pd.to_datetime(data_['date_max'], format='%Y-%m-%d')

    return data_

def append_variables(path='',file='',data=None,category=''):

    data_ = data

    text = open(path+'/'+file).read()
    text = re.sub('\\pm','~',text)
    text = re.sub(r'[\\£@#${}^]','',text)

    #Parse top line for variables
    variables = text.split("\n")[0].split(' & ')

    for i in range(2,len(variables),1):
        values = []
        errors = []
        for line in text.split('\n'):

            if line.split('-runNumber')[0] == category:

                value = line.split(' & ')[i]

                if '--' in value: value = '0'

                if len(value.split('~')) == 2:
                    values.append(float(value.split('~')[0]))
                    errors.append(float(value.split('~')[1]))
                elif len(value.split('~')) == 1:
                    values.append(float(value))
                    errors.append(0.0)

        if 'mc' in variables[i]:
            variables[i] = variables[i].replace("mc","MC")
        if '/' in variables[i]:
            variables[i] = variables[i].replace("/","_")

        data_[variables[i]] = values
        data_[variables[i]+'_err'] = errors

    return data_


def plot_stability( xData = None, datavalues = None, mcvalues = None, dataerrors = None, mcerrors = None, label = '', category = '', path = "", evenX = False, xVar = ''):

    left, width = 0.1, 1.0
    bottom, height = 0.1, 0.5
    rect_hist = [left+width+0.01, bottom, 0.2, height]
    rect_plot = [left, bottom, width, height]

    nullfmt = NullFormatter()

    fig = plt.figure()

    ax_plot = plt.axes(rect_plot)
    ax_hist = plt.axes(rect_hist)

    ax_hist.yaxis.set_major_formatter(nullfmt)

    xPlaceholder = range(1,1+len(xData),1)
    if evenX:
        ax_plot.errorbar(xPlaceholder,datavalues,yerr=dataerrors,capthick=0,marker='o',ms=5,ls='None',)
    else:
        ax_plot.errorbar(xData,datavalues,yerr=dataerrors,capthick=0,marker='o',ms=5,ls='None',)

    # customise the axes
    xDataVar = xData.name
    if xDataVar == 'time':
        ax_plot.xaxis.set_minor_locator(dates.DayLocator(interval=2))
        ax_plot.xaxis.set_minor_formatter(dates.DateFormatter('%d\n%a'))
        ax_plot.xaxis.set_major_locator(dates.MonthLocator())
        ax_plot.xaxis.set_major_formatter(dates.DateFormatter('\n\n\n%b\n%Y'))
    elif (xDataVar == 'run_max' or xDataVar == 'run_min') and not evenX:
        majorLocator = MultipleLocator(125)
        minorLocator = MultipleLocator(62.5)
        ax_plot.xaxis.set_major_locator(majorLocator)
        ax_plot.xaxis.set_minor_locator(minorLocator)
        majorFormatter = FormatStrFormatter('%d')
        ax_plot.xaxis.set_major_formatter(majorFormatter)
        xlabels = ax_plot.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=10)
    elif (xDataVar == 'run_max' or xDataVar == 'run_min') and evenX:
        majorLocator = MultipleLocator(2)
        minorLocator = MultipleLocator(1)
        ax_plot.xaxis.set_major_locator(majorLocator)
        ax_plot.xaxis.set_minor_locator(minorLocator)
        xlabels = ax_plot.get_xticks().tolist()

        for i in range(2,len(xData),2):
            xlabels[i/2+1] = xData.tolist()[i-1]

        for i in range(len(xlabels)):
            if xlabels[i] < 200000: xlabels[i] = ''

        ax_plot.set_xticklabels(xlabels)
        xlabels = ax_plot.get_xticklabels()
        plt.setp(xlabels, rotation=90, fontsize=10)

    ax_plot.xaxis.grid(True, which="minor")
    ax_plot.yaxis.grid()
    ax_hist.xaxis.grid(True, which="minor")
    ax_hist.yaxis.grid()
    ax_hist.xaxis.set_ticks([])

    #Get and set the limits for the histogram
    ymin = round(datavalues.min()) - 1
    ymax = round(datavalues.max()) + 1
    ax_plot.set_ylim((ymin,ymax))
    ax_hist.set_ylim((ymin,ymax))

    ax_plot.set_ylabel(label)

    nbin = 20
    y,_,_ = ax_hist.hist(datavalues, bins=nbin, orientation='horizontal', histtype='stepfilled', alpha=0.6)
    hmax = y.max()
    ax_hist.set_xlim((0,hmax*1.1))

    ax_plot.set_title(region_labels[category] + "    " + label)

    #Annotate with mean and std dev
    npVals = np.asarray(datavalues)
    ax_hist.annotate('Mean = {:3.3f}'.format(np.mean(npVals)),(hmax/6,ymin-(ymax-ymin)*0.1),fontsize=11,annotation_clip=False,xycoords='data')
    ax_hist.annotate('Std dev. = {:3.3f}'.format(np.std(npVals)),(hmax/6,ymin-(ymax-ymin)*0.175),fontsize=11,annotation_clip=False,xycoords='data')

    #Add line for the MC
    if (len(mcvalues) > 0):
        xmin,xmax = 0,1
        ax_plot.errorbar(xData,mcvalues,yerr=mcerrors,capthick=0,marker='o',ms=5,ls='None',c='Red')

        xNP = np.asarray(xData.tolist())
        mcNP = np.asarray(mcvalues.tolist())
        mcErrNP = np.asarray(mcerrors.tolist())

        ax_plot.fill_between(xNP,mcNP-mcErrNP,mcNP+mcErrNP,alpha=0.3,edgecolor='red', facecolor='red')


#        ax_plot.axhline(y=mcvalue+mcerror,
#                        xmin=xmin,xmax=xmax,c='red',linewidth=1,clip_on=True,linestyle='dashed')
#        ax_plot.axhline(y=mcvalue-mcerror,
#                        xmin=xmin,xmax=xmax,c='red',linewidth=1,clip_on=True,linestyle='dashed')
#        ax_plot.axhline(y=mcvalue,
#                        xmin=xmin,xmax=xmax,c='blue',linewidth=1,clip_on=True,linestyle='solid',label='MC')
#        ax_hist.annotate('MC = {:3.3f} $\pm$ {:3.3f}'.format(mcvalue,mcerror),(hmax/6,ymin-(ymax-ymin)*0.25),fontsize=11,annotation_clip=False,xycoords='data')

    #Legend
    # legend = ax_plot.legend(loc='lower center',numpoints=1)
    # if (mcvalue > -999):
    #     legend.get_texts()[1].set_text('Data')
    # else:
    #     legend.get_texts()[0].set_text('Data')

    #Save
    if evenX:
        xDataVar = xDataVar + '_even'

    if not os.path.exists(path+'/'+xDataVar):
        os.makedirs(path+'/'+xDataVar)

    pdfName = category + '_' + re.sub(r'[\\@#$/{}^]','', label) + '.pdf'
    pdfName = re.sub(r'[ ]','_',pdfName)

    print 'Saving plot: ' + path + xDataVar + '/' + pdfName
    for fr in format_fig_output :
        plt.savefig( path + xDataVar + '/' + pdfName,
                     format='pdf',orientation='landscape',
                     dpi=900,papertype='a0',pad_inches=0.1,
                     bbox_inches='tight')

    plt.close(fig)
