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

doShifts = True


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

# TODO reorder to only have one in memory at once
counters = {}
shiftTypes = ['lep','trig','pu','fake','ElectronEn','MuonEn','TauEn','JetEn','UnclusteredEn']
shiftTypes = ['trig','pu','fake','ElectronEn','MuonEn','TauEn','JetEn','UnclusteredEn'] # lep broken somehow
shifts = ['']
for s in shiftTypes:
    shifts += [s+'Up',s+'Down']
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

def getCount(counters,sig,directory,shift=''):
    tot, totErr = counters[sig+shift].getCount(sig,directory)
    return (tot,totErr)

def getBackgroundCount(counters,directory,datadriven=False,shift=''):
    tot = 0
    totErr2 = 0
    if datadriven:
        dirdirs = directory.split('/')
        for s in samples:
            val,err = getCount(counters,s,'/'.join(['4P0F']+dirdirs),shift=shift)
            tot += val
            totErr2 += err**2
        for s in samples+['data']:
            for reg in ['3P1F','2P2F','1P3F','0P4F']:
                val,err = getCount(counters,s,'/'.join([reg]+dirdirs),shift=shift)
                tot += val
                totErr2 += err**2
    else:
        for s in allsamples:
            val,err = getCount(counters,s,directory,shift=shift)
            tot += val
            totErr2 += err**2
    return (tot,totErr2**0.5)

def getAlphaCount(counters,directory,datadriven=False,alphaOnly=False,shift=''):
    mc_side       = getBackgroundCount(counters,'new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_mw         = getBackgroundCount(counters,'new/massWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    #mc_all        = getBackgroundCount(counters,'new/allMassWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alpha         = divWithError(mc_mw,mc_side)
    if alphaOnly: return alpha[0]
    data_allside  = getCount(counters,'data','new/allSideband/{0}'.format(directory))
    data_exp      = prodWithError(data_allside,alpha)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    return (abs(data_exp[0]),abs(data_allside[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha

def getAlphaPrimeCount(counters,directory,datadriven=False,alphaOnly=False,shift=''):
    mc_side       = getBackgroundCount(counters,'new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_allSide    = getBackgroundCount(counters,'new/allSideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alpha         = divWithError(mc_allSide,mc_side)
    if alphaOnly: return alpha[0]
    data_side     = getCount(counters,'data','new/sideband/{0}'.format(directory))
    data_exp      = prodWithError(data_side,alpha)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    return (abs(data_exp[0]),abs(data_side[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha


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

def getRecoChans(mode):
    # find out what reco/gen channels can exist for this mode
    recoChans = set()
    for gen in genRecoMap:
        if len(gen)!=4: continue # only 4l allowed here
        s = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
        if not s: continue
        recoChans.update(genRecoMap[gen])
    return recoChans

# get shifted first
allShifts = {}
for shift in shifts[1:]:
    signals = ['HppHmm{0}GeV'.format(mass) for mass in masses]
    # counters
    counters = {}
    for s in allsamples:
        counters[s+shift] = Counter('Hpp4l')
        counters[s+shift].addProcess(s,sigMap[s],shift=shift)

    for s in signals:
        counters[s+shift] = Counter('Hpp4l')
        counters[s+shift].addProcess(s,sigMap[s],signal=True,shift=shift)

    counters['data'+shift] = Counter('Hpp4l')
    counters['data'+shift].addProcess('data',sigMap['data'],shift=shift)
    # get values
    allShifts[shift] = {}
    for mode in modes:
        allShifts[shift][mode] = {}
        for mass in masses:
            logging.info('Getting shift {0}: {1} - {2} GeV'.format(shift,mode,mass))
            recoChans = getRecoChans(mode)
            signals = ['HppHmm{0}GeV'.format(mass)]
            results = {}
            for reco in recoChans:
                recoSB = reco+'_SB'
                results[reco] = {}
                results[recoSB] = {}
                # for 100%, get num taus, for benchmarks, based on reco
                hpphmm = 'hpp{0}hmm{1}'.format(modeMap[mode][0],modeMap[mode][1])
                # get the shifts effect on alpha
                results[reco]['alpha'] = getAlphaCount(counters,'{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=2,shift=shift,alphaOnly=True)
                results[recoSB]['alpha'] = getAlphaPrimeCount(counters,'{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=2,shift=shift,alphaOnly=True)
                # PP
                for proc in signals:
                    results[reco][proc] = 0.
                    results[recoSB][proc] = 0.
                    for gen in genRecoMap:
                        if len(gen)!=4: continue # 3 for AP, 4 for PP
                        if reco not in genRecoMap[gen]: continue
                        scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                        results[reco][proc] += scale*getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen),shift=shift)[0]
                        results[recoSB][proc] += scale*getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen),shift=shift)[0]
            allShifts[shift][mode][mass] = results


signals = ['HppHmm{0}GeV'.format(mass) for mass in masses]
# counters
counters = {}
for s in allsamples:
    counters[s] = Counter('Hpp4l')
    counters[s].addProcess(s,sigMap[s])

for s in signals:
    counters[s] = Counter('Hpp4l')
    counters[s].addProcess(s,sigMap[s],signal=True)

counters['data'] = Counter('Hpp4l')
counters['data'].addProcess('data',sigMap['data'])
for mode in modes:
    for mass in masses:
        logging.info('Producing datacard for {0} - {1} GeV'.format(mode,mass))
        results = {}
        limits = Limits()
    
        limits.addEra('13TeV80X')
        limits.addAnalysis('Hpp4l')
        
        recoChans = getRecoChans(mode)
        for reco in recoChans:
            limits.addChannel(reco)
            limits.addChannel(reco+'_SB')

        signals = ['HppHmm{0}GeV'.format(mass)]
        for sig in signals:
            limits.addProcess(sig,signal=True)
        
        for background in backgrounds:
            limits.addProcess(background)

        # set values and stat error
        staterr = {}
        uncerr = {x:{} for x in shiftTypes}
        era = '13TeV80X'
        analysis = 'Hpp4l'
        for reco in recoChans:
            recoSB = reco+'_SB'
            results[reco] = {}
            results[recoSB] = {}
            # for 100%, get num taus, for benchmarks, based on reco
            hpphmm = 'hpp{0}hmm{1}'.format(modeMap[mode][0],modeMap[mode][1])
            if len(backgrounds)==1 and backgrounds[0] == 'datadriven':
                value,side,alpha,err = getAlphaCount(counters,'{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=3)
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
                for unc in shiftTypes:
                    err = (abs(allShifts[unc+'Up'][mode][mass][reco]['alpha']-alpha)+abs(allShifts[unc+'Down'][mode][mass][reco]['alpha']-alpha))/2.
                    if err and alpha: uncerr[unc][(('datadriven',),(era,),(analysis,),(reco,))] = 1+err/alpha
                # sideband values
                value,side,alpha,err = getAlphaPrimeCount(counters,'{0}/{1}/{2}'.format(mass,hpphmm,reco),datadriven=reco.count('t')>=3)
                limits.setExpected('datadriven',era,analysis,recoSB,value)
                limits.addSystematic('alphaSB_{era}_{analysis}_{channel}'.format(era=era,analysis=analysis,channel=recoSB),
                                     'gmN {0}'.format(int(side)),
                                     systematics={(('datadriven',),(era,),(analysis,),(recoSB,)):alpha})
                if value: staterr[(('datadriven',),(era,),(analysis,),(recoSB,))] = 1+err/value
                results[recoSB]['alpha'] = alpha
                results[recoSB]['alphaError'] = err
                results[recoSB]['expected'] = value
                results[recoSB]['sideCount'] = side
                # get the shifts effect on alpha
                for unc in shiftTypes:
                    err = (abs(allShifts[unc+'Up'][mode][mass][recoSB]['alpha']-alpha)+abs(allShifts[unc+'Down'][mode][mass][recoSB]['alpha']-alpha))/2.
                    if err and alpha: uncerr[unc][(('datadriven',),(era,),(analysis,),(recoSB,))] = 1+err/alpha
            else:
                for proc in backgrounds:
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphmm,reco))
                    limits.setExpected(proc,era,analysis,reco,value)
                    if value: staterr[((proc,),(era,),(analysis,),(reco,))] = 1+err/value
            for proc in signals:
                totalValue = 0.
                err2 = 0.
                totalValueSB = 0.
                err2SB = 0.
                for gen in genRecoMap:
                    if len(gen)!=4: continue # only 4l allowed here
                    if reco not in genRecoMap[gen]: continue
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen))
                    scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                    totalValue += scale*value
                    err2 += (scale*err)**2
                    # sideband
                    value,err = getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphmm,reco,gen))
                    totalValueSB += scale*value
                    err2SB += (scale*err)**2
                limits.setExpected(proc,era,analysis,reco,totalValue)
                if totalValue: staterr[((proc,),(era,),(analysis,),(reco,))] = 1.+err2**0.5/totalValue
                results[reco]['pp'] = totalValue
                results[reco]['ppError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-totalValue)+abs(allShifts[unc+'Down'][mode][mass][reco][proc]-totalValue))/2.
                    if err and totalValue: uncerr[unc][((proc,),(era,),(analysis,),(reco,))] = 1+err/totalValue
                # sideband
                limits.setExpected(proc,era,analysis,recoSB,totalValueSB)
                if totalValueSB: staterr[((proc,),(era,),(analysis,),(recoSB,))] = 1.+err2SB**0.5/totalValueSB
                results[recoSB]['pp'] = totalValueSB
                results[recoSB]['ppError'] = err2SB**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-totalValueSB)+abs(allShifts[unc+'Down'][mode][mass][recoSB][proc]-totalValueSB))/2.
                    if err and totalValueSB: uncerr[unc][((proc,),(era,),(analysis,),(recoSB,))] = 1+err/totalValueSB
            obs = getCount(counters,'data','new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphmm,reco))
            limits.setObserved(era,analysis,reco,obs)
            results[reco]['observed'] = obs[0]
            # sideband
            obs = getCount(counters,'data','new/allSideband/{0}/{1}/{2}'.format(mass,hpphmm,recoSB))
            limits.setObserved(era,analysis,recoSB,obs)
            results[recoSB]['observed'] = obs[0]
            dumpResults(results,'Hpp4l','{0}/{1}'.format(mode,mass))

        # systematics
        addUncertainties(limits,staterr,uncerr,recoChans,signals,backgrounds,4)

        # print the datacard
        directory = 'datacards/{0}/{1}'.format('Hpp4l',mode)
        python_mkdir(directory)
        limits.printCard('{0}/{1}.txt'.format(directory,mass))
