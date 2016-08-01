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
from DevTools.Limits.higgsUncertainties import addUncertainties

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

doShifts = False


# define cards to create
modes = ['ee100','em100','et100','mm100','mt100','tt100','BP1','BP2','BP3','BP4']
masses = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]

cats = getCategories('Hpp4l')
catLabels = getCategoryLabels('Hpp4l')
subCatChannels = getSubCategories('Hpp4l')
subCatLabels = getSubCategoryLabels('Hpp4l')
chans = getChannels('Hpp4l')
chanLabels = getChannelLabels('Hpp4l')
genRecoMap = getGenRecoChannelMap('Hpp4l')
sigMap = getSigMap('Hpp4l')
sigMapDD = getSigMap('Hpp4l',datadriven=True)

scales = {}
for mode in modes:
    scales[mode] = getScales(mode)

samples = ['TTV','VVV','ZZ']
allsamples = ['TT','TTV','Z','WZ','VVV','ZZ']
signals = ['HppHmm{0}GeV'.format(mass) for mass in masses]
backgrounds = ['datadriven']

datadrivenSamples = []
for s in samples + ['data']:
    datadrivenSamples += sigMap[s]

counters = {}
shifts = ['','lepUp','lepDown','trigUp','trigDown','puUp','puDown','fakeUp','fakeDown']
shiftTypes = ['lep','trig','pu','fake']
if not doShifts:
    shifts = ['']
    shiftTypes = []
for shift in shifts:
    for s in allsamples:
        counters[s+shift] = Counter('Hpp4l')
        counters[s+shift].addProcess(s,sigMap[s],shift=shift)
    
    for s in signals:
        counters[s+shift] = Counter('Hpp4l')
        counters[s+shift].addProcess(s,sigMap[s],signal=True,shift=shift)
    
    counters['data'+shift] = Counter('Hpp4l')
    counters['data'+shift].addProcess('data',sigMap['data'],shift=shift)

def getCount(sig,directory,shift=''):
    tot, totErr = counters[sig+shift].getCount(sig,directory)
    return (tot,totErr)

def getBackgroundCount(directory,datadriven=False,shift=''):
    tot = 0
    totErr2 = 0
    if datadriven:
        dirdirs = directory.split('/')
        for s in samples:
            val,err = getCount(s,'/'.join(['4P0F']+dirdirs),shift=shift)
            tot += val
            totErr2 += err**2
        for s in samples+['data']:
            for reg in ['3P1F','2P2F','1P3F','0P4F']:
                val,err = getCount(s,'/'.join([reg]+dirdirs),shift=shift)
                tot += val
                totErr2 += err**2
    else:
        for s in allsamples:
            val,err = getCount(s,directory,shift=shift)
            tot += val
            totErr2 += err**2
    return (tot,totErr2**0.5)

def getAlphaCount(directory,datadriven=False,alphaOnly=False,shift=''):
    mc_side       = getBackgroundCount('new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_mw         = getBackgroundCount('new/massWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    #mc_all        = getBackgroundCount('new/allMassWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alpha         = divWithError(mc_mw,mc_side)
    if alphaOnly: return alpha[0]
    data_allside  = getCount('data','new/allSideband/{0}'.format(directory))
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
        results = {}
        limits = Limits()
    
        limits.addEra('13TeV80X')
        limits.addAnalysis('Hpp4l')
        
        # find out what reco/gen channels can exist for this mode
        recoChans = set()
        for gen in genRecoMap:
            if len(gen)!=4: continue # only 4l allowed here
            s = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
            if not s: continue
            recoChans.update(genRecoMap[gen])
        for reco in recoChans: limits.addChannel(reco)

        signals = ['HppHmm{0}GeV'.format(mass)]
        for sig in signals:
            limits.addProcess(sig,signal=True)
        
        for background in backgrounds:
            limits.addProcess(background)

        # set values and stat error
        staterr = {}
        uncerr = {x:{} for x in shiftTypes}
        for era in ['13TeV80X']:
            for analysis in ['Hpp4l']:
                for reco in recoChans:
                    results[reco] = {}
                    # for 100%, get num taus, for benchmarks, based on reco
                    hpphmm = 'hpp{0}hmm{1}'.format(modeMap[mode][0],modeMap[mode][1])
                    if len(backgrounds)==1 and backgrounds[0] == 'datadriven':
                        value,side,alpha,err = getAlphaCount('{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=3)
                        limits.setExpected('datadriven',era,analysis,reco,value)
                        limits.addSystematic('alpha_{era}_{analysis}_{channel}'.format(era=era,analysis=analysis,channel=reco),
                                             'gmN {0}'.format(int(side)),
                                             systematics={(('datadriven',),(era,),(analysis,),(reco,)):alpha})
                        if value: staterr[(('datadriven',),(era,),(analysis,),(reco,))] = 1+err/value
                        results[reco]['alpha'] = alpha
                        results[reco]['alphaError'] = err
                        results[reco]['expected'] = value
                        results[reco]['sideCount'] = side
                        # get the shifts effect on alpha
                        allAlpha = {}
                        for shift in shifts[1:]:
                            allAlpha[shift] = getAlphaCount('{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=2,shift=shift,alphaOnly=True)
                        for unc in shiftTypes:
                            err = abs((allAlpha[unc+'Up']+allAlpha[unc+'Down'])/2.-alpha)
                            if err and alpha: uncerr[unc][(('datadriven',),(era,),(analysis,),(reco,))] = 1+err/alpha
                    else:
                        for proc in backgrounds:
                            value,err = getCount(proc,'new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphmm,reco))
                            limits.setExpected(proc,era,analysis,reco,value)
                            if value: staterr[((proc,),(era,),(analysis,),(reco,))] = 1+err/value
                    for proc in signals:
                        totalValue = 0.
                        err2 = 0.
                        allShifts = {x:0. for x in shifts[1:]}
                        for gen in genRecoMap:
                            if len(gen)!=4: continue # only 4l allowed here
                            if reco not in genRecoMap[gen]: continue
                            value,err = getCount(proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen))
                            scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                            totalValue += scale*value
                            err2 += (scale*err)**2
                            for shift in shifts[1:]:
                                allShifts[shift] += scale*getCount(proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen),shift=shift)[0]
                        limits.setExpected(proc,era,analysis,reco,totalValue)
                        if totalValue: staterr[((proc,),(era,),(analysis,),(reco,))] = 1.+err2**0.5/totalValue
                        results[reco]['pp'] = totalValue
                        results[reco]['ppError'] = err2**0.5
                        # get the shifts effect on signal
                        for unc in shiftTypes:
                            err = abs((allShifts[unc+'Up']+allShifts[unc+'Down'])/2.-totalValue)
                            if err and totalValue: uncerr[unc][((proc,),(era,),(analysis,),(reco,))] = 1+err/totalValue
                    obs = getCount('data','new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphmm,reco))
                    limits.setObserved(era,analysis,reco,obs)
                    results[reco]['observed'] = obs[0]
                    dumpResults(results,'Hpp4l','{0}/{1}'.format(mode,mass))

        # systematics
        addUncertainties(limits,staterr,uncerr,recoChans,signals,backgrounds,4)

        # print the datacard
        directory = 'datacards/{0}/{1}'.format('Hpp4l',mode)
        python_mkdir(directory)
        limits.printCard('{0}/{1}.txt'.format(directory,mass))
