import os
import sys
import logging

import ROOT

from DevTools.Limits.Limits import Limits
from DevTools.Plotter.NtupleWrapper import NtupleWrapper
from DevTools.Utilities.utilities import *
from DevTools.Plotter.threePhotonUtils import *
import DevTools.Limits.Models as Models

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

# control
blind = True
selection = 'pass_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90 && gg12_mass>100 && g1_passPreselection>0 && g2_passPreselection>0 && g1_mvaNonTrigValues>0. && g2_mvaNonTrigValues>0. && g3_mvaNonTrigValues>0.'
scalefactor = '*'.join(['genWeight','pileupWeight',])
binned = True
addSignal = True
wsname = 'w'
doParametric = False
binning = [10,0,500]

# setup
sampleMap = getSampleMap()

backgrounds = ['TT','TTG','VVG','Z','G','GG','QCD']
data = ['data']
#signals = ['HToAG_250_1','HToAG_250_30','HToAG_250_150']
signals = ['HToAG_250_150']

wrappers = {}
for proc in backgrounds+signals+data:
    for sample in sampleMap[proc]:
        wrappers[sample] = NtupleWrapper('ThreePhoton',sample,new=True,version='80X')

def getBinned(proc,**kwargs):
    scalefactor = kwargs.pop('scalefactor','*'.join(['genWeight','pileupWeight',]))
    hists = ROOT.TList()
    for sample in sampleMap[proc]:
        hist = wrappers[sample].getTempHist2D(sample,selection,'1' if proc=='data' else scalefactor,'gg13_mass','gg23_mass',binning,binning)
        hists.Add(hist)
    if hists.IsEmpty():
        hist = 0
    else:
        hist = hists[0].Clone('h_{0}'.format(proc))
        hist.Reset()
        hist.Merge(hists)
    return hist

def getUnbinned(proc):
    return ROOT.RooDataSet()

def sumHists(name,*hists):
    histlist = ROOT.TList()
    for hist in hists:
        histlist.Add(hist)
    hist = histlist[0].Clone(name)
    hist.Reset()
    hist.Merge(histlist)
    return hist

# load histograms
histMap = {}
for proc in backgrounds+signals:
    logging.info('Getting {0}'.format(proc))
    histMap[proc] = getBinned(proc)
logging.info('Getting observed')
if blind:
    samples = backgrounds
    if addSignal: samples = backgrounds + signals
    hists = []
    for proc in samples:
        hists += [histMap[proc]]
    hist = sumHists('obs',*hists)
    for b in range(hist.GetNbinsX()+1):
        val = int(hist.GetBinContent(b))
        if val<0: val = 0
        err = val**0.5
        hist.SetBinContent(b,val)
        #hist.SetBinError(b,err)
    histMap['data'] = hist
else:
    hist = getBinned('data')
    histMap['data'] = hist

# create limit object
limits = Limits(wsname)

limits.addEra('Run2016')
limits.addAnalysis('ThreePhoton')
limits.addChannel('ggg')

if doParametric:
    limits.addMH(*binning[1:])
    limits.addX(*binning[1:])

if doParametric:
    limits.addProcess('sig',signal=True)
    limits.addProcess('bg')
else:
    for signal in signals:
        limits.addProcess(signal,signal=True)
    for background in backgrounds:
        limits.addProcess(background)

era = 'Run2016'
analysis = 'ThreePhoton'
reco = 'ggg'
if doParametric:
    # no parametric for 2d fit yet
    pass

else:
    for proc in backgrounds:
        limits.setExpected(proc,era,analysis,reco,histMap[proc])
    for proc in signals:
        limits.setExpected(proc,era,analysis,reco,histMap[proc])

limits.setObserved(era,analysis,reco,histMap['data'])

# add systematics
systproc = tuple([proc for proc in signals + backgrounds if 'datadriven' not in proc])

############
### Lumi ###
############
# lumi 2.3% for 2015 and 2.5% for 2016
# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
logging.info('Adding lumi systematic')
lumisyst = {
    (systproc,(era,),('all',),('all',)): 1.025,
}
limits.addSystematic('lumi','lnN',systematics=lumisyst)

##############
### Pileup ###
##############
puMapUp = {}
puMapDown = {}
for proc in backgrounds+signals:
    logging.info('Getting {0} PU up'.format(proc))
    puMapUp[proc] = getBinned(proc,scalefactor='*'.join(['genWeight','pileupWeightUp',]))
    logging.info('Getting {0} PU down'.format(proc))
    puMapDown[proc] = getBinned(proc,scalefactor='*'.join(['genWeight','pileupWeightDown',]))
pusyst = {}
for proc in systproc:
    pusyst[((proc,),(era,),(analysis,),(reco,))] = (puMapUp[proc],puMapDown[proc])
limits.addSystematic('pu','shape',systematics=pusyst)

# print the datacard
directory = 'datacards_shape/{0}'.format('ThreePhoton')
python_mkdir(directory)
datacard = '{0}/ggg_mod.txt'.format(directory)
processes = ['sig','bg'] if doParametric else signals+backgrounds
limits.printCard(datacard,processes=processes,blind=False,saveWorkspace=doParametric)

