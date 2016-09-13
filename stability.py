#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import sys, getopt
from stability  import ParseTable as pt
from shutil     import copyfile
from optparse   import OptionParser
from pprint     import pprint
import pandas   as pd

def get_options():
    parser = OptionParser()
    parser.add_option("-r", "--run-range",
                      dest="inRunRange", default='',
                      help="""
                      Input configuration file that contains the run ranges.
                      This file is produced by the Zfitter (please check the documentation)
                      and can be found in data/runRange/ directory
                      """,
                      metavar="FILE")

    parser.add_option("-s", "--stability",
                      dest="inStability", default='stability.tex',
                      help="""
                      stability tex file produced by the ZFitter
                      """,
                      metavar="FILE")

    parser.add_option("--iovs",
                      dest="iovs_file", default='',
                      help="""
                      List of IOV's to be drawn on the history plots
                      The file should be structured as folowing :
                      | run | date | time | info |
                      """,
                      metavar="FILE")

    parser.add_option("-i", "--inv-mass",
                      dest="invMass", default='invMassSC',
                      help="""
                      Invariant mass variable.
                      """)

    parser.add_option('-x','--xVar',
                      dest='xVar',default='',
                      help='''
                      Specify whether the plot should be over run range, or regions
                      ''')

    parser.add_option('-d','--outputName',
                      dest='outputName',default='',
                      help="""
                      Output directory for the plots
                      """)

    return parser.parse_args()


if __name__ == '__main__':

    (opt, args) =  get_options()

    assert opt.inStability != '', 'No stability file input given, use -s <stability.tex>'
    assert opt.inRunRange  != '' or opt.xVar != '', 'No runranges or variable to scan over specified'
    assert opt.invMass     != '', 'No invariant mass name specified, use -i <invariant mass>'
    assert opt.outputName      != '', 'No output name specified'

    if opt.inRunRange != '':
        print "Run range file to be used is "   ,opt.inRunRange
    else:
        print "Scan variable is ",opt.xVar
    print "Stability file to be used is "   ,opt.inStability
    print "Invariant mass that was used is ",opt.invMass
    print "Plot output directory is plot-stability/",opt.outputName

    if opt.xVar != '':
        print "x-axis variable is ", opt.xVar
    else:
        print "x-axis variables are the times/runNumbers"

    if not os.path.exists('plot-stability/'):
        os.makedirs('plot-stability/')
    if not os.path.exists('data-stability/'):
        os.makedirs('data-stability/')

    #Rename and copy the datafiles into data/
    runRangeFile  = opt.inRunRange.split('/')[-1]
    stabilityFile = opt.inRunRange.split('/')[-1].replace('.dat','')+'_'+opt.invMass+'_'+opt.xVar+'_stability.tex'

    data_path = 'data-stability/'
    plot_path = 'plot-stability/' + runRangeFile.replace('.dat','') + '/' + opt.invMass + '/'

    if opt.inRunRange != data_path + runRangeFile:
        copyfile(opt.inRunRange , data_path +'/'+ runRangeFile)
    if opt.inStability != data_path + stabilityFile:
        copyfile(opt.inStability, data_path +'/'+ stabilityFile)
        print "Stability .tex file is now at " + data_path +'/'+ stabilityFile

    if opt.inRunRange != '':
        runRangeFile  = opt.inRunRange.split('/')[-1]
        if opt.inRunRange != data_path + runRangeFile:
            copyfile(opt.inRunRange , data_path +'/'+ runRangeFile)
            print "Runrange file is now at " + data_path +'/'+ runRangeFile


    #Reading regions from the table file
    regions = pt.read_regions_from_table(path=data_path,tableFile=stabilityFile,xVar=opt.xVar)
    print "categories :: ", regions

    # reading iov's
    iovs = None
    if os.path.exists(opt.iovs_file):
        iovs =  pd.read_csv( opt.iovs_file ,sep=' ', names = ['run', 'date', 'time', 'playload', 'info'])
    print "Starting plotmaking..."
    for region in regions:
        print "Category: ",region

        if opt.xVar != '':
            d=pt.parse_table_over_regions(path=data_path, tableFile=stabilityFile,category=region,xVar=opt.xVar)
        else:
            #Get runrange and time info from the the runranges file
            d = pt.read_run_range(path=data_path,file=runRangeFile)
            #Get variables information from the stability monitoring .tex file
            d = pt.append_variables(path=data_path,file=stabilityFile,data=d,category=region)

        #Get variables to make plots of (data, not mc or err vars)
        variables = []
        if opt.xVar != '':
            xVars = [ opt.xVar + '_min', opt.xVar + '_max' ]
        else:
            xVars = [ 'Nevents'     ,
                      'UnixTime'    ,
                      'run_number'  ,
                      'UnixTime_min',
                      'UnixTime_max',
                      'run_min'     ,
                      'run_max'     ,
                      'date_min'    ,
                      'date_max'    ,
                      'time'        ]
        for label in d.columns.values.tolist():
            if "MC" not in label and label not in xVars and "_err" not in label:
                variables.append(label)

        #Loop over the vars
        for var in variables:
            #Get associated monte carlo info, or a placeholder
            varmc = var.replace("data","MC")
            if 'MC' not in varmc:
                print "[WARNING] MC counterpart not found for ", var
                mcvalues = []
                mcerrors = []
            else:
                mcvalues = d[varmc]
                mcerrors = d[varmc+'_err']


            if opt.xVar == '':
                #Switches on whether the datapoints are evenly distributed along x
                evenXs = [False,True]
                #Plot as function of date or run numbers
                timevars = ['run_min']
                for timevar in timevars:
                    for evenX in evenXs:
                        pt.plot_stability( xData = d[timevar], datavalues = d[var],
                                           dataerrors = d[var+'_err'], mcvalues = mcvalues,
                                           mcerrors = mcerrors, label = pt.var_labels[var],
                                           category = region, path=plot_path, evenX = evenX,
                                           xVar=opt.xVar, iovs=iovs, var=var)
            else:
                if opt.xVar in var : continue
                pt.plot_stability( xData  = d[opt.xVar], xData_err=d[opt.xVar + '_err'], datavalues = d[var],
                                   dataerrors = d[var+'_err'], mcvalues = mcvalues,
                                   mcerrors   = mcerrors, label = pt.var_labels[var],
                                   category   = region  , path=plot_path, evenX = False,
                                   xVar       = opt.xVar, iovs=iovs, var=var)
