import os
import sys
import logging
import itertools

import ROOT
ROOT.gROOT.SetBatch()

from DevTools.Limits.Limits import Limits
from DevTools.Plotter.NtupleWrapper import NtupleWrapper
from DevTools.Utilities.utilities import *
from DevTools.Plotter.haaUtils import *
import DevTools.Limits.Models as Models

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

###############
### Control ###
###############
blind = True
selection = '1'
scalefactor = '*'.join(['genWeight','pileupWeight','triggerEfficiency'])
addSignal = False
wsname = 'w'
mmbinning = [290,1,30]
mtbinning = [600,0,60]
hbinning  = [100,0,1000]

#############
### Setup ###
#############
sampleMap = getSampleMap()

#backgrounds = ['JPsi','Upsilon', 'W', 'Z', 'TT', 'WW', 'WZ', 'ZZ']
#backgrounds = ['W', 'Z', 'TT', 'WW', 'WZ', 'ZZ']
#backgrounds = ['W', 'Z', 'TT']
backgrounds = ['datadriven']
data = ['data']
signame = 'HToAAH{h}A{a}'

hmasses = [125,300,750]
amasses = [5,7,9,11,13,15,17,19,21]

signals = [signame.format(h=h,a=a) for h in hmasses for a in amasses]
#signals = [signame.format(h=125,a=15)]

wrappers = {}
for proc in backgrounds+signals+data:
    if proc=='datadriven': continue
    for sample in sampleMap[proc]:
        wrappers[sample] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X')

#################
### Utilities ###
#################
def getBinned(procname,**kwargs):
    scalefactor = kwargs.pop('scalefactor','1' if procname=='data' else '*'.join(['genWeight','pileupWeight','triggerEfficiency']))

    proc = procname

    # use region D for datadriven background and scale by 0.5
    if procname=='datadriven':
        selection = '(am1_isolation>0.15 || am2_isolation>0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT<0.5'
        proc = 'data'
        scalefactor = '0.5'
    else:
        selection = 'am1_isolation<0.15 && am2_isolation<0.15 && ath_byVLooseIsolationMVArun2v1DBoldDMwLT>0.5'

    hists = ROOT.TList()
    for sample in sampleMap[proc]:
        hist = wrappers[sample].getTempHist(sample,selection,scalefactor,'amm_mass',mmbinning)
        hists.Add(hist)
    if hists.IsEmpty():
        hist = 0
    else:
        hist = hists[0].Clone('h_{0}'.format(procname))
        hist.Reset()
        hist.Merge(hists)
    return hist

def sumHists(name,*hists):
    histlist = ROOT.TList()
    for hist in hists:
        histlist.Add(hist)
    hist = histlist[0].Clone(name)
    hist.Reset()
    hist.Merge(histlist)
    return hist

##############################
### Create/read histograms ###
##############################

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

#####################
### Create Limits ###
#####################
limits = Limits(wsname)

limits.addEra('Run2016')
limits.addAnalysis('HAA')
limits.addChannel('mmmt')

era = 'Run2016'
analysis = 'HAA'
reco = 'mmmt'

for signal in signals:
    limits.addProcess(signal,signal=True)
for background in backgrounds:
    limits.addProcess(background)

for proc in backgrounds:
    limits.setExpected(proc,era,analysis,reco,histMap[proc])
for proc in signals:
    limits.setExpected(proc,era,analysis,reco,histMap[proc])

limits.setObserved(era,analysis,reco,histMap['data'])

#########################
### Add uncertainties ###
#########################


systproc = tuple([proc for proc in signals + backgrounds if 'datadriven' not in proc])

############
### stat ###
############

def getStat(hist,direction):
    newhist = hist.Clone('{0}{1}'.format(hist.GetName(),direction))
    nb = hist.GetNbinsX()*hist.GetNbinsY()
    for b in range(nb+1):
        val = hist.GetBinContent(b+1)
        err = hist.GetBinError(b+1)
        newval = val+err if direction=='Up' else val-err
        if newval<0: newval = 0
        newhist.SetBinContent(b+1,newval)
        newhist.SetBinError(b+1,0)
    return newhist

logging.info('Adding stat systematic')
statMapUp = {}
statMapDown = {}
for proc in systproc:
    statMapUp[proc] = getStat(histMap[proc],'Up')
    statMapDown[proc] = getStat(histMap[proc],'Down')

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
logging.info('Adding pileup systematic')
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

######################
### Print datacard ###
######################
directory = 'datacards_shape/{0}'.format('MuMuTauTau')
python_mkdir(directory)
datacard = '{0}/mmmt_1d_{1}'.format(directory,'binned')
processes = {}
for signal in signals:
    processes[signal] = [signal]+backgrounds
limits.printCard(datacard,processes=processes,blind=False)
