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
regions = {
    'A': '(am1_isolation<0.15 && am2_isolation<0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT>0.5',
    'B': '(am1_isolation<0.15 && am2_isolation<0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT<0.5',
    'C': '(am1_isolation>0.15 || am2_isolation>0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT>0.5',
    'D': '(am1_isolation>0.15 || am2_isolation>0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT<0.5',
}

doParametric = False
blind = True
selection = regions['A']
scalefactor = '*'.join(['genWeight','pileupWeight','triggerEfficiency'])
addSignal = True
signalParams = {'h': 125, 'a': 15}
wsname = 'w'
mmbinning = [260,4,30]
mtbinning = [600,0,60]
hbinning  = [100,0,1000]

#############
### Setup ###
#############
sampleMap = getSampleMap()

#backgrounds = ['JPsi','Upsilon', 'W', 'Z', 'TT', 'WW', 'WZ', 'ZZ']
#backgrounds = ['W', 'Z', 'TT', 'WW', 'WZ', 'ZZ']
#backgrounds = ['W', 'Z', 'TT']
backgrounds = ['datadriven','TT']
data = ['data']
signame = 'HToAAH{h}A{a}'

hmasses = [125,300,750]
hmasses = [125]
amasses = [5,7,9,11,13,15,17,19,21]
amasses = [5,11,15,21]

signals = [signame.format(h=h,a=a) for h in hmasses for a in amasses]
#signals = [signame.format(h=125,a=15)]
signalToAdd = signame.format(**signalParams)

wrappers = {}
for proc in backgrounds+signals+data:
    if proc=='datadriven': continue
    for sample in sampleMap[proc]:
        wrappers[sample] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X')

#################
### Utilities ###
#################
def getBinned(proc,**kwargs):
    sf = kwargs.pop('scalefactor','1' if proc=='data' else scalefactor)
    sel = kwargs.pop('selection',selection)

    hists = ROOT.TList()
    for sample in sampleMap[proc]:
        hist = wrappers[sample].getTempHist(sample,sel,sf,'amm_mass',mmbinning)
        hists.Add(hist)
    if hists.IsEmpty():
        hist = 0
    else:
        hist = hists[0].Clone('h_{0}'.format(proc))
        hist.Reset()
        hist.Merge(hists)
    return hist

def getDatadriven(**kwargs):
    sf = kwargs.pop('scalefactor',scalefactor)

    # region D
    # scale down by 0.5 for now
    thisSel = regions['D']

    hists = ROOT.TList()

    #first get data
    data = getBinned('data',selection=thisSel,scalefactor='0.5',**kwargs)
    hists.Add(data)

    # get all MC backgrounds and subtract
    for proc in backgrounds:
        if 'datadriven' in proc: continue
        hist = getBinned(proc,selection=thisSel,scalefactor='-0.5*{0}'.format(sf),**kwargs)
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

def getSpline(histMap,h):
    # initial fit
    results = {}
    results[h] = {}
    for a in amasses:
        ws = ROOT.RooWorkspace('sig')
        ws.factory('x[{0}, {1}]'.format(*mmbinning[1:]))
        model = Models.Voigtian('sig',
            mean  = [a,0,30],
            width = [0.01*a,0,5],
            sigma = [0.01*a,0,5],
        )
        model.build(ws, 'sig')
        hist = histMap[signame.format(h=h,a=a)]
        results[h][a] = model.fit(ws, hist, '{0}_{1}'.format(h,a), save=True)

    # create model
    for a in amasses:
        print h, a, results[h][a]
    model = Models.VoigtianSpline(splinename.format(h=h),
        **{
            'masses' : amasses,
            'means'  : [results[h][a]['mean_{0}_{1}'.format(h,a)] for a in amasses],
            'widths' : [results[h][a]['width_{0}_{1}'.format(h,a)] for a in amasses],
            'sigmas' : [results[h][a]['sigma_{0}_{1}'.format(h,a)] for a in amasses],
        }
    )
    integrals = [histMap[signame.format(h=h,a=a)].Integral() for a in amasses]
    integral = np.mean(integrals)
    model.setIntegral(integral)

    return model

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

if doParametric:

    limits.addMH(*mmbinning[1:])
    limits.addX(*mmbinning[1:])
    for h in hmasses:
        limits.addProcess(splinename.format(h=h),signal=True)
    for background in backgrounds:
        limits.addProcess(background)
    
    # add models
    for h in hmasses:
        model = getSpline(histMap,h)
        limits.setExpected(splinename.format(h=h),era,analysis,reco,model)
    
    # add histograms
    for bg in backgrounds:
        limits.setExpected(bg,era,analysis,reco,histMap[bg])
    
    # get roodatahist
    limits.setObserved(era,analysis,reco,histMap['data'])

else:

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
allproc = tuple([proc for proc in signals + backgrounds])
systsplineproc = tuple([proc for proc in signalSplines + backgrounds if 'datadriven' not in proc])
allsplineproc = tuple([proc for proc in signalSplines + backgrounds])
bgproc = tuple([proc for proc in backgrounds])
sigsplineproc = tuple([proc for proc in signalSplines])
sigproc = tuple([proc for proc in signals])


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
for proc in backgrounds+signals:
    statMapUp[proc] = getStat(histMap[proc],'Up')
    statMapDown[proc] = getStat(histMap[proc],'Down')
statsyst = {}
for proc in bgproc:
    statsyst[((proc,),(era,),(analysis,),(reco,))] = (statMapUp[proc],statMapDown[proc])
if doParametric:
    for h in hmasses:
        statsyst[((splinename.format(h=h),),(era,),(analysis,),(reco,))] = (getSpline(statMapUp,h),getSpline(statMapDown,h))
else:
    for proc in sigproc:
        statsyst[((proc,),(era,),(analysis,),(reco,))] = (statMapUp[proc],statMapDown[proc])
limits.addSystematic('stat_{process}_{channel}','shape',systematics=statsyst)

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
    if proc=='datadriven':
        puMapUp[proc] = getDatadriven(scalefactor='*'.join(['genWeight','pileupWeightUp','triggerEfficiency']))
    else:
        puMapUp[proc] = getBinned(proc,scalefactor='*'.join(['genWeight','pileupWeightUp','triggerEfficiency']))
    logging.info('Getting {0} PU down'.format(proc))
    if proc=='datadriven':
        puMapDown[proc] = getDatadriven(scalefactor='*'.join(['genWeight','pileupWeightDown','triggerEfficiency']))
    else:
        puMapDown[proc] = getBinned(proc,scalefactor='*'.join(['genWeight','pileupWeightDown','triggerEfficiency']))
pusyst = {}
for proc in bgproc:
    pusyst[((proc,),(era,),(analysis,),(reco,))] = (puMapUp[proc],puMapDown[proc])
if doParametric:
    for h in hmasses:
        pusyst[((splinename.format(h=h),),(era,),(analysis,),(reco,))] = (getSpline(puMapUp,h),getSpline(puMapDown,h))
else:
    for proc in sigproc:
        pusyst[((proc,),(era,),(analysis,),(reco,))] = (puMapUp[proc],puMapDown[proc])
limits.addSystematic('pu','shape',systematics=pusyst)

######################
### Print datacard ###
######################
directory = 'datacards_shape/{0}'.format('MuMuTauTau')
python_mkdir(directory)
datacard = '{0}/mmmt_1d_{1}'.format(directory,'parametric' if doParametric else 'binned')
processes = {}
if doParametric:
    for h in hmasses:
        processes[signame.format(h=h,a='X')] = [splinename.format(h=h)] + backgrounds
else:
    for signal in signals:
        processes[signal] = [signal]+backgrounds
limits.printCard(datacard,processes=processes,blind=False,saveWorkspace=doParametric)
