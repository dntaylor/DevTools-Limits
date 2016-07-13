import os
import sys
import logging
import ROOT
import numpy as np
import math

from DevTools.Limits.Limits import Limits
from DevTools.Utilities.utilities import *
from DevTools.Plotter.Counter import Counter
from DevTools.Plotter.higgsUtilities import *

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


# define cards to create
modes = ['ee100','em100','et100','mm100','mt100','tt100','BP1','BP2','BP3','BP4']
masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]

cats = getCategories('Hpp3l')
catLabels = getCategoryLabels('Hpp3l')
subCatChannels = getSubCategories('Hpp3l')
subCatLabels = getSubCategoryLabels('Hpp3l')
chans = getChannels('Hpp3l')
chanLabels = getChannelLabels('Hpp3l')
genRecoMap = getGenRecoChannelMap('Hpp3l')
sigMap = getSigMap('Hpp3l')
sigMapDD = getSigMap('Hpp3l',datadriven=True)

scales = {}
for mode in modes:
    scales[mode] = getScales(mode)

samples = ['TTV','VVV','ZZ']
allsamples = ['TT','TTV','Z','WZ','VVV','ZZ']
signals = ['HppHm{0}GeV'.format(mass) for mass in masses]
backgrounds = ['datadriven']

datadrivenSamples = []
for s in samples + ['data']:
    datadrivenSamples += sigMap[s]

counters = {}
for s in allsamples:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s])

for s in signals:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s],signal=True)

counters['data'] = Counter('Hpp3l')
counters['data'].addProcess('data',sigMap['data'])

def getCount(sig,directory):
    tot, totErr = counters[sig].getCount(sig,directory)
    return (tot,totErr)

def getBackgroundCount(directory):
    tot = 0
    totErr2 = 0
    for s in allsamples:
        sname = s.replace('all','')
        val,err = getCount(sname,directory)
        tot += val
        totErr2 += err**2
    return (tot,totErr2**0.5)

def getAlphaCount(directory):
    mc_side       = getBackgroundCount('old/sideband/{0}'.format(directory))
    mc_mw         = getBackgroundCount('old/massWindow/{0}'.format(directory))
    #mc_all        = getBackgroundCount('old/allMassWindow/{0}'.format(directory))
    data_allside  = getCount('data','old/allSideband/{0}'.format(directory))
    alpha         = divWithError(mc_mw,mc_side)
    data_exp      = prodWithError(data_allside,alpha)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    return (abs(data_exp[0]),abs(data_allside[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha

# TODO, think if this is what we want
modeMap = {
    'ee100': [0,0],
    'em100': [0,0],
    'et100': [1,1],
    'mm100': [0,0],
    'mt100': [1,1],
    'tt100': [2,2],
    'BP1'  : [2,2],
    'BP2'  : [2,2],
    'BP3'  : [2,2],
    'BP4'  : [2,2],
}

for mode in modes:
    for mass in masses:
        logging.info('Producing datacard for {0} - {1} GeV'.format(mode,mass))
        limits = Limits()
    
        limits.addEra('13TeV')
        limits.addAnalysis('Hpp3l')
        
        # find out what reco/gen channels can exist for this mode
        recoChans = set()
        for gen in genRecoMap:
            if len(gen)!=3: continue # only 3l allowed here
            s = scales[mode].scale_Hpp3l(gen[:2],gen[2:])
            if not s: continue
            recoChans.update(genRecoMap[gen])
        for reco in recoChans: limits.addChannel(reco)

        signals = ['HppHm{0}GeV'.format(mass)]
        for sig in signals:
            limits.addProcess(sig,signal=True)
        
        for background in backgrounds:
            limits.addProcess(background)

        # set values and stat error
        staterr = {}
        for era in ['13TeV']:
            for analysis in ['Hpp3l']:
                for reco in recoChans:
                    # for 100%, get num taus, for benchmarks, based on reco
                    hpphm = 'hpp{0}'.format(modeMap[mode][0])
                    if len(backgrounds)==1 and backgrounds[0] == 'datadriven':
                        value,side,alpha,err = getAlphaCount('{0}/{1}/{2}'.format(mass,hpphm,reco))
                        limits.setExpected('datadriven',era,analysis,reco,value)
                        limits.addSystematic('alpha_{era}_{analysis}_{channel}'.format(era=era,analysis=analysis,channel=reco),
                                             'gmN {0}'.format(int(side)),
                                             systematics={(('datadriven',),(era,),(analysis,),(reco,)):alpha})
                        if value: staterr[(('datadriven',),(era,),(analysis,),(reco,))] = 1+err/value
                    else:
                        for proc in backgrounds:
                            value,err = getCount(proc,'old/allMassWindow/{0}/{1}/{2}'.format(mass,hpphm,reco))
                            limits.setExpected(proc,era,analysis,reco,value)
                            if value: staterr[((proc,),(era,),(analysis,),(reco,))] = 1+err/value
                    for proc in signals:
                        totalValue = 0.
                        err2 = 0.
                        for gen in genRecoMap:
                            if len(gen)!=3: continue # only 3l allowed here
                            if reco not in genRecoMap[gen]: continue
                            value,err = getCount(proc,'old/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                            scale = scales[mode].scale_Hpp3l(gen[:2],gen[2:])
                            totalValue += scale*value
                            err2 += (scale*err)**2
                        limits.setExpected(proc,era,analysis,reco,totalValue)
                        if totalValue: staterr[((proc,),(era,),(analysis,),(reco,))] = 1.+err2**0.5/totalValue
                    obs = getCount('data','old/allMassWindow/{0}/{1}/{2}'.format(mass,hpphm,reco))
                    limits.setObserved(era,analysis,reco,obs)

        # systematics
        # stat errs
        limits.addSystematic('stat_{process}','lnN',systematics=staterr)

        systproc = tuple([proc for proc in signals + backgrounds if proc!='datadriven'])

        # lumi 2.7% for 2015 and 2016
        lumisyst = {
            (systproc,('13TeV',),('all',),('all',)): 1.027,
        }
        limits.addSystematic('lumi_{era}','lnN',systematics=lumisyst)

        # electron id 2%/leg
        elecsyst = {}
        for c in range(3):
            systChans = tuple([chan for chan in recoChans if chan.count('e')==c+1])
            if not systChans: continue
            elecsyst[(systproc,('13TeV',),('Hpp3l',),systChans)] = 1.+math.sqrt((c+1)*0.02**2)
        if elecsyst: limits.addSystematic('elec_id','lnN',systematics=elecsyst)

        # muon id 1+0.5%/leg
        muonsyst = {}
        for c in range(3):
            systChans = tuple([chan for chan in recoChans if chan.count('m')==c+1])
            if not systChans: continue
            muonsyst[(systproc,('13TeV',),('Hpp3l',),systChans)] = 1.+math.sqrt((c+1)*(0.01**2 + 0.005**2))
        if muonsyst: limits.addSystematic('muon_id','lnN',systematics=muonsyst)

        # taus id 6%
        tausyst = {}
        for c in range(3):
            systChans = tuple([chan for chan in recoChans if chan.count('t')==c+1])
            if not systChans: continue
            tausyst[(systproc,('13TeV',),('Hpp3l',),systChans)] = 1.+math.sqrt((c+1)*(0.06**2))
        if tausyst: limits.addSystematic('tau_id','lnN',systematics=tausyst)

        # print the datacard
        directory = 'datacards/{0}/{1}'.format('Hpp3l',mode)
        python_mkdir(directory)
        limits.printCard('{0}/{1}.txt'.format(directory,mass))
