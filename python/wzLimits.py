import os
import sys
import logging
import ROOT
import math

from DevTools.Limits.Limits import Limits
from DevTools.Plotter.Counter import Counter
from DevTools.Utilities.utilities import *
from DevTools.Plotter.wzUtilities import sigMap


chans = ['eee','eem','mme','mmm']
samples = ['TTV','ZG','VVV','ZZ','WZ']
allsamples = ['W','TT','Z','WW','TTV','VVV','ZZall','WZall']

def getCount(counters,sig,directory,shift=''):
    tot, totErr = counters[sig+shift].getCount(sig,directory)
    return (tot,totErr)

def getDatadrivenCount(counters,directory,shift=''):
    tot = 0
    totErr2 = 0
    dirdirs = directory.split('/')
    for s in samples+['data']:
        for reg in ['2P1F','1P2F','0P3F']:
            val,err = getCount(counters,s,'/'.join([reg]+dirdirs),shift=shift)
            tot += val
            totErr2 += err**2
    return (tot,totErr2**0.5)

counters = {}
for s in samples+['data']:
    counters[s] = Counter('WZ')
    counters[s].addProcess(s,sigMap[s])

counts = {}
for chan in chans:
    counts[chan] = {}
    for sample in samples+['data']:
        counts[chan][sample] = getCount(counters,sample,'3P0F/default/{0}'.format(chan))
    counts[chan]['datadriven'] = getDatadrivenCount(counters,'default/{0}'.format(chan))

era = '13TeV80X'
analysis = 'WZ'

limits = Limits()
limits.addEra(era)
limits.addAnalysis(analysis)
for chan in chans:
    limits.addChannel(chan)
for sig in ['WZ']:
    limits.addProcess(sig,signal=True)
for bg in ['datadriven','ZZ','TTV','VVV','ZG']:
    limits.addProcess(bg)

staterr = {}
for chan in chans:
    for process in samples+['datadriven']:
        val,err = counts[chan][process]
        limits.setExpected(process,era,analysis,chan,val)
        if val: staterr[((process,),(era,),(analysis,),(chan,))] = 1+err/val
    limits.setObserved('13TeV80X','WZ',chan,counts[chan]['data'][0])

limits.addSystematic('stat_{process}_{channel}','lnN',systematics=staterr)

limits.printCard('wz.txt',blind=False)
for chan in chans:
    limits.printCard('wz_{0}.txt'.format(chan),channels=[chan],blind=False)
