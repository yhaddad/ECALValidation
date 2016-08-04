#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import sys, getopt, time
from monitoring  import ParseTable as pt
from shutil      import copyfile
from optparse    import OptionParser
from pprint      import pprint

def get_options():
    parser = OptionParser()
    parser.add_option("-t", "--table",
                      dest="intable", default='stability.tex',
                      help="""
                      table containing the scan results from ZFitter
                      """,
                      metavar="FILE")

    parser.add_option("-v", "--variable",
                      dest="variable", default='absEta',
                      help="""
                      Variable scanned by ZFitter
                      """)
    parser.add_option("-l", "--variable-title",
                      dest="vtitle", default='|\eta|',
                      help="""
                      Latex type title. to be intruduced without '$ $'
                      example : -l |\eta|
                      """)
    parser.add_option("-d", "--make-date-dir",
                      action="store_true",dest="make-date-dir",default=True,
                      help="""
                      Store the results and the tables in a dated directories
                      """)

    return parser.parse_args()


if __name__ == '__main__':
    (opt, args) =  get_options()
    assert opt.intable   != '', 'No scan table given , please use -t <scan_table.tex>'
    assert opt.variable  != '', 'No variable specified , please use -v <variable>'

    print "-- Scan table       : " ,opt.intable
    print "-- Scanned variable : " ,opt.variable
    print "-- Variable title   : " ,opt.vtitle

    outdirname = 'plot-scan/' + time.strftime("%d-%m-%Y") + '/'
    # if not os.path.exists(outdirname):
        # os.makedirs(outdirname)

    #Rename and copy the datafiles into data/
    data_path  = opt.intable.split('/')[-1]
    intable    = opt.intable.split('/')[-1]
    if opt.intable  != outdirname + intable:
        copyfile(opt.intable , intable +'/'+ intable)
    if opt.inStability != data_path + stabilityFile:
        copyfile(opt.inStability, data_path +'/'+ stabilityFile)

    print "They can be found in data/ as " + runRangeFile + " and " + stabilityFile
    """
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)

    print "Starting plotmaking..."
    print "categories :: "

    regions = pt.read_regions_from_table(path=data_path,tableFile=stabilityFile)
    print "Regions are:"
    for region in regions: print region

    if overRegions:
        print "Plots are over regions"
    else:
        print "Plots are over run ranges"

    for category in regions:
        print "Beginning category ", category,
        if "gold" in category:
            print " ...skipping: gold"
            continue
        if "bad"  in category:
            print " ...skipping: bad"
            continue
        print

        #Get runrange and time info from the the runranges file
        d = pt.read_run_range(path=data_path,file=runRangeFile)

        #Get variables information from the stability monitoring .tex file
        d = pt.append_variables(path=data_path,file=stabilityFile,data=d,category=category)

        #Get variables to make plots of (data, not mc or err vars)
        variables = []
        timeVars = ['Nevents'     ,
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
            if "MC" not in label and label not in timeVars and "_err" not in label:
                variables.append(label)

        #Loop over the vars
        for var in variables:
            #Get associated monte carlo info, or a placeholder
            varmc = var.replace("data","MC")
            if 'MC' not in varmc:
                print "[WARNING] MC counterpart not found for ", var
                mcvalue = -999
                mcerror = -999
            else:
                mcvalue = d[varmc][0]
                mcerror = d[varmc+'_err'][0]

            #Switches on whether the datapoints are evenly distributed along x
            evenXs = [False,True]

            #Plot as function of date or run numbers
            timevars = ['run_min','run_max','time']
            for timevar in timevars:
                for evenX in evenXs:
                    pt.plot_stability( time = d[timevar], datavalues = d[var],
                                       dataerrors = d[var+'_err'], mcvalue = mcvalue,
                                       mcerror = mcerror, label = pt.var_labels[var],
                                       category = category, path=plot_path, evenX = evenX )
    """
