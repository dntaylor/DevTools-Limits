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

blind = False
doShifts = True
readUncerr = False # read from file rather than compute on the fly

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

samples = ['TTV','VVV','ZZ','WZ']
allsamples = ['W','T','TT','TTV','Z','WW','VVV','ZZ','WZ']
signalsAP = ['HppHm{0}GeV'.format(mass) for mass in masses]
signalsPP = ['HppHmm{0}GeV'.format(mass) for mass in masses]
signalsPPR = ['HppHmmR{0}GeV'.format(mass) for mass in masses]
backgrounds = ['datadriven']

datadrivenSamples = []
for s in samples + ['data']:
    datadrivenSamples += sigMap[s]

counters = {}
#shiftTypes = ['lep','trig','pu','fake','btag','ElectronEn','MuonEn','TauEn','JetEn','UnclusteredEn']
shiftTypes = ['lep','trig','pu','fake','btag','ElectronEn','MuonEn','TauEn','JetEn']
#shiftTypes = ['lep','trig','pu','fake','btag','ElectronEn','MuonEn','TauEn']
shifts = ['']
for s in shiftTypes:
    shifts += [s+'Up',s+'Down']
if not doShifts:
    shifts = ['']
    shiftTypes = []

def getCount(counters,sig,directory,shift=''):
    tot, totErr = counters[sig+shift].getCount(sig,directory)
    return (tot,totErr)

def getBackgroundCount(counters,directory,datadriven=False,shift=''):
    tot = 0
    totErr2 = 0
    if datadriven:
        dirdirs = directory.split('/')
        for s in samples:
            val,err = getCount(counters,s,'/'.join(['3P0F']+dirdirs),shift=shift)
            tot += val
            totErr2 += err**2
        for s in samples+['data']:
            for reg in ['2P1F','1P2F','0P3F']:
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
    #mc_side       = getBackgroundCount(counters,'new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    #mc_mw         = getBackgroundCount(counters,'new/massWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_mw          = getBackgroundCount(counters,'new/massWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_allmw       = getBackgroundCount(counters,'new/allMassWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    #mc_all        = getBackgroundCount(counters,'new/allMassWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alpha         = divWithError(mc_allmw,mc_mw)
    if abs(alpha[0]) < abs(alpha[1]): alpha = (alpha[1], alpha[1])
    if alphaOnly: return abs(alpha[0])
    #data_allside  = getCount(counters,'data','new/allSideband/{0}'.format(directory))
    data_mw  = getCount(counters,'data','new/massWindow/{0}'.format(directory))
    data_exp      = prodWithError(data_mw,alpha)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    #return (abs(data_exp[0]),abs(data_allside[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha
    return (abs(data_exp[0]),abs(data_mw[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha

def getAlphaPrimeCount(counters,directory,datadriven=False,alphaOnly=False,shift=''):
    mc_side       = getBackgroundCount(counters,'new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_allSide    = getBackgroundCount(counters,'new/allSideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alpha         = divWithError(mc_allSide,mc_side)
    if abs(alpha[0]) < abs(alpha[1]): alpha = (alpha[1], alpha[1])
    if alphaOnly: return abs(alpha[0])
    data_side     = getCount(counters,'data','new/sideband/{0}'.format(directory))
    data_exp      = prodWithError(data_side,alpha)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    return (abs(data_exp[0]),abs(data_side[0]),abs(alpha[0]),abs(alpha[1])) # fix for negative alpha

def getDualAlphaCount(counters,directory,datadriven=False,alphaOnly=False,shift=''):
    mc_side       = getBackgroundCount(counters,'new/sideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_allSide    = getBackgroundCount(counters,'new/allSideband/{0}'.format(directory),datadriven=datadriven,shift=shift)
    mc_allmw      = getBackgroundCount(counters,'new/allMassWindow/{0}'.format(directory),datadriven=datadriven,shift=shift)
    alphaSR       = divWithError(mc_allmw,mc_side)
    alphaSB       = divWithError(mc_allSide,mc_side)
    if abs(alphaSR[0]) < abs(alphaSR[1]): alphaSR = (alphaSR[1], alphaSR[1])
    if abs(alphaSB[0]) < abs(alphaSB[1]): alphaSB = (alphaSB[1], alphaSB[1])
    if alphaOnly: return (abs(alphaSR[0]),abs(alphaSB[0]))
    data_side     = getCount(counters,'data','new/sideband/{0}'.format(directory))
    dataSR_exp    = prodWithError(data_side,alphaSR)
    dataSB_exp    = prodWithError(data_side,alphaSB)
    # return data_exp, data_sideband, alpha, alpha stat uncertainty
    return (abs(dataSR_exp[0]),abs(dataSB_exp[0]),abs(data_side[0]),abs(alphaSR[0]),abs(alphaSR[1]),abs(alphaSB[0]),abs(alphaSB[1])) # fix for negative alpha

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
        if len(gen)==3:
            s = scales[mode].scale_Hpp3l(gen[:2],gen[2:])
        else:
            s = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
        if not s: continue
        recoChans.update(genRecoMap[gen])
    return recoChans

# get shifted first
allShifts = {}
for shift in shifts[1:]:
    if readUncerr: continue
    signalsAP = ['HppHm{0}GeV'.format(mass) for mass in masses]
    signalsPP = ['HppHmm{0}GeV'.format(mass) for mass in masses]
    signalsPPR = ['HppHmmR{0}GeV'.format(mass) for mass in masses]
    # counters
    counters = {}
    for s in allsamples:
        counters[s+shift] = Counter('Hpp3l')
        counters[s+shift].addProcess(s,sigMap[s],shift=shift)
    
    for s in signalsAP:
        counters[s+shift] = Counter('Hpp3l')
        counters[s+shift].addProcess(s,sigMap[s],signal=True,shift=shift)
    
    for s in signalsPP:
        counters[s+shift] = Counter('Hpp3l')
        counters[s+shift].addProcess(s,sigMap[s],signal=True,shift=shift)

    for s in signalsPPR:
        counters[s+shift] = Counter('Hpp3l')
        counters[s+shift].addProcess(s,sigMap[s],signal=True,shift=shift)

    counters['data'+shift] = Counter('Hpp3l')
    counters['data'+shift].addProcess('data',sigMap['data'],shift=shift)
    # get values
    allShifts[shift] = {}
    for mode in modes:
        allShifts[shift][mode] = {}
        for mass in masses:
            logging.info('Getting shift {0}: {1} - {2} GeV'.format(shift,mode,mass))
            recoChans = getRecoChans(mode)
            signalsAP = ['HppHm{0}GeV'.format(mass)]
            signalsPP = ['HppHmm{0}GeV'.format(mass)]
            signalsPPR = ['HppHmmR{0}GeV'.format(mass)]
            results = {}
            for reco in recoChans:
                recoSB = reco + '_SB'
                results[reco] = {}
                results[recoSB] = {}
                # for 100%, get num taus, for benchmarks, based on reco
                hpphm = 'hpp{0}'.format(modeMap[mode][0])
                # get the shifts effect on alpha
                results[reco]['alpha'], results[recoSB]['alpha'] = getDualAlphaCount(counters,'{0}/{1}/{2}'.format(mass,hpphm,reco),datadriven=True,shift=shift,alphaOnly=True)#reco.count('t')>=2 or reco[-1]=='t',shift=shift,alphaOnly=True)
                # AP
                for proc in signalsAP:
                    results[reco][proc] = 0.
                    results[recoSB][proc] = 0.
                    for gen in genRecoMap:
                        if len(gen)!=3: continue # 3 for AP, 4 for PP
                        if reco not in genRecoMap[gen]: continue
                        scale = scales[mode].scale_Hpp3l(gen[:2],gen[2:])
                        results[reco][proc] += scale*getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
                        results[recoSB][proc] += scale*getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
                # PP
                for proc in signalsPP:
                    results[reco][proc] = 0.
                    results[recoSB][proc] = 0.
                    for gen in genRecoMap:
                        if len(gen)!=4: continue # 3 for AP, 4 for PP
                        if reco not in genRecoMap[gen]: continue
                        scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                        results[reco][proc] += scale*getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
                        results[recoSB][proc] += scale*getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
                # PPR
                for proc in signalsPPR:
                    results[reco][proc] = 0.
                    results[recoSB][proc] = 0.
                    for gen in genRecoMap:
                        if len(gen)!=4: continue # 3 for AP, 4 for PP
                        if reco not in genRecoMap[gen]: continue
                        scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                        results[reco][proc] += scale*getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
                        results[recoSB][proc] += scale*getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen),shift=shift)[0]
            allShifts[shift][mode][mass] = results

signalsAP = ['HppHm{0}GeV'.format(mass) for mass in masses]
signalsPP = ['HppHmm{0}GeV'.format(mass) for mass in masses]
signalsPPR = ['HppHmmR{0}GeV'.format(mass) for mass in masses]
# counters
counters = {}
for s in allsamples:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s])

for s in signalsAP:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s],signal=True)

for s in signalsPP:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s],signal=True)

for s in signalsPPR:
    counters[s] = Counter('Hpp3l')
    counters[s].addProcess(s,sigMap[s],signal=True)

counters['data'] = Counter('Hpp3l')
counters['data'].addProcess('data',sigMap['data'])
for mode in modes:
    for mass in masses:
        logging.info('Producing datacard for {0} - {1} GeV'.format(mode,mass))
        results = {}
        limits = Limits()
    
        limits.addEra('Era13TeV2016')
        limits.addAnalysis('Hpp3l')
        limits.addAnalysis('Hpp3lAP')
        limits.addAnalysis('Hpp3lPP')
        limits.addAnalysis('Hpp3lPPR')
        
        recoChans = getRecoChans(mode)
        for reco in recoChans:
            limits.addChannel(reco)
            limits.addChannel(reco+'_SB')

        signalsAP = ['HppHm{0}GeV'.format(mass)]
        signalsPP = ['HppHmm{0}GeV'.format(mass)]
        signalsPPR = ['HppHmmR{0}GeV'.format(mass)]
        for sig in signalsAP+signalsPP+signalsPPR:
            limits.addProcess(sig,signal=True)
        
        for background in backgrounds:
            limits.addProcess(background)

        # set values and stat error
        staterr = {}
        uncerr = {x:{} for x in shiftTypes}
        uncerr_store = {x:{} for x in shiftTypes}
        era = 'Era13TeV2016'
        analysis = 'Hpp3l'

        # read uncerrtainties from file
        if readUncerr:
            uncerr_store = readResults(analysis,'shiftUncertainties_{0}_{1}'.format(mode,mass))

        for reco in recoChans:
            recoSB = reco + '_SB'
            results[reco] = {}
            results[recoSB] = {}
            # for 100%, get num taus, for benchmarks, based on reco
            hpphm = 'hpp{0}'.format(modeMap[mode][0])
            if len(backgrounds)==1 and backgrounds[0] == 'datadriven':
                # dual alpha
                valueSR,valueSB,side,alphaSR,errSR,alphaSB,errSB = getDualAlphaCount(counters,'{0}/{1}/{2}'.format(mass,hpphm,reco),datadriven=True)#reco.count('t')>=2 or reco[-1]=='t')
                limits.setExpected('datadriven',era,analysis,reco,valueSR)
                limits.setExpected('datadriven',era,analysis+'AP',reco,valueSR)
                limits.setExpected('datadriven',era,analysis+'PP',reco,valueSR)
                limits.setExpected('datadriven',era,analysis+'PPR',reco,valueSR)
                limits.setExpected('datadriven',era,analysis,recoSB,valueSB)
                limits.setExpected('datadriven',era,analysis+'AP',recoSB,valueSB)
                limits.setExpected('datadriven',era,analysis+'PP',recoSB,valueSB)
                limits.setExpected('datadriven',era,analysis+'PPR',recoSB,valueSB)
                limits.addSystematic('alpha_{era}_{analysis}_{channel}'.format(era=era,analysis=analysis,channel=reco),
                                     'gmN {0}'.format(int(side)),
                                     systematics={
                                         (('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(reco,)):alphaSR,
                                         (('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(recoSB,)):alphaSB,
                                     })
                if alphaSR: staterr[(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(reco,))] = 1+errSR/alphaSR
                if alphaSB: staterr[(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(recoSB,))] = 1+errSB/alphaSB
                results[reco]['alpha'] = alphaSR
                results[reco]['alphaError'] = errSR
                results[reco]['expected'] = valueSR
                results[reco]['sideCount'] = side
                results[recoSB]['alpha'] = alphaSB
                results[recoSB]['alphaError'] = errSB
                results[recoSB]['expected'] = valueSB
                results[recoSB]['sideCount'] = side
                # get the shifts effect on alpha
                for unc in shiftTypes:
                    if readUncerr:
                        if 'datadriven_{0}'.format(reco) in uncerr_store[unc]:
                            uncerr[unc][(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(reco,))] = 1+uncerr_store[unc]['datadriven_{0}'.format(reco)]
                        if 'datadriven_{0}'.format(recoSB) in uncerr_store[unc]:
                            uncerr[unc][(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(recoSB,))] = 1+uncerr_store[unc]['datadriven_{0}'.format(recoSB)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][reco]['alpha']-alphaSR)+abs(allShifts[unc+'Down'][mode][mass][reco]['alpha']-alphaSR))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][reco]['alpha']-allShifts[unc+'Down'][mode][mass][reco]['alpha']))/2.
                        if err and alphaSR:
                            uncerr[unc][(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(reco,))] = min([1+err/alphaSR,2])
                            if not readUncerr: uncerr_store[unc]['datadriven_{0}'.format(reco)] = err/alphaSR
                        #err = (abs(allShifts[unc+'Up'][mode][mass][recoSB]['alpha']-alphaSB)+abs(allShifts[unc+'Down'][mode][mass][recoSB]['alpha']-alphaSB))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][recoSB]['alpha']-allShifts[unc+'Down'][mode][mass][recoSB]['alpha']))/2.
                        if err and alphaSB:
                            uncerr[unc][(('datadriven',),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(recoSB,))] = min([1+err/alphaSB,2])
                            uncerr_store[unc]['datadriven_{0}'.format(recoSB)] = err/alphaSB

            else:
                for proc in backgrounds:
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphm,reco))
                    limits.setExpected(proc,era,analysis,reco,value)
                    limits.setExpected(proc,era,analysis+'AP',reco,value)
                    limits.setExpected(proc,era,analysis+'PP',reco,value)
                    limits.setExpected(proc,era,analysis+'PPR',reco,value)
                    if value: staterr[((proc,),(era,),(analysis,analysis+'AP',analysis+'PP',analysis+'PPR',),(reco,))] = 1+err/value
            # AP
            for proc in signalsAP:
                totalValue = 0.
                err2 = 0.
                totalValueSB = 0.
                err2SB = 0.
                for gen in genRecoMap:
                    if len(gen)!=3: continue # 3 for AP, 4 for PP
                    if reco not in genRecoMap[gen]: continue
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    scale = scales[mode].scale_Hpp3l(gen[:2],gen[2:])
                    totalValue += scale*value
                    err2 += (scale*err)**2
                    # sb values
                    value,err = getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    totalValueSB += scale*value
                    err2SB += (scale*err)**2
                limits.setExpected(proc,era,analysis,reco,totalValue)
                limits.setExpected(proc,era,analysis+'AP',reco,totalValue)
                if totalValue: staterr[((proc,),(era,),(analysis,analysis+'AP',),(reco,))] = 1.+err2**0.5/totalValue
                results[reco]['ap'] = totalValue
                results[reco]['apError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,reco) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'AP',),(reco,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,reco)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-totalValue)+abs(allShifts[unc+'Down'][mode][mass][reco][proc]-totalValue))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-allShifts[unc+'Down'][mode][mass][reco][proc]))/2.
                        if err and totalValue:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'AP',),(reco,))] = min([1+err/totalValue,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,reco)] = err/totalValue
                # sideband components
                limits.setExpected(proc,era,analysis,recoSB,totalValueSB)
                limits.setExpected(proc,era,analysis+'AP',recoSB,totalValueSB)
                if totalValueSB: staterr[((proc,),(era,),(analysis,analysis+'AP',),(recoSB,))] = 1.+err2SB**0.5/totalValueSB
                results[recoSB]['ap'] = totalValueSB
                results[recoSB]['apError'] = err2SB**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,recoSB) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'AP',),(recoSB,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,recoSB)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-totalValueSB)+abs(allShifts[unc+'Down'][mode][mass][recoSB][proc]-totalValueSB))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-allShifts[unc+'Down'][mode][mass][recoSB][proc]))/2.
                        if err and totalValueSB:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'AP',),(recoSB,))] = min([1+err/totalValueSB,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,recoSB)] = err/totalValueSB
            # PP
            for proc in signalsPP:
                totalValue = 0.
                err2 = 0.
                totalValueSB = 0.
                err2SB = 0.
                for gen in genRecoMap:
                    if len(gen)!=4: continue # 3 for AP, 4 for PP
                    if reco not in genRecoMap[gen]: continue
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                    totalValue += scale*value
                    err2 += (scale*err)**2
                    # sb values
                    value,err = getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    totalValueSB += scale*value
                    err2SB += (scale*err)**2
                limits.setExpected(proc,era,analysis,reco,totalValue)
                limits.setExpected(proc,era,analysis+'PP',reco,totalValue)
                if totalValue: staterr[((proc,),(era,),(analysis,analysis+'PP',),(reco,))] = 1.+err2**0.5/totalValue
                results[reco]['pp'] = totalValue
                results[reco]['ppError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,reco) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'PP',),(reco,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,reco)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-totalValue)+abs(allShifts[unc+'Down'][mode][mass][reco][proc]-totalValue))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-allShifts[unc+'Down'][mode][mass][reco][proc]))/2.
                        if err and totalValue:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'PP',),(reco,))] = min([1+err/totalValue,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,reco)] = err/totalValue
                # sideband components
                limits.setExpected(proc,era,analysis,recoSB,totalValueSB)
                limits.setExpected(proc,era,analysis+'PP',recoSB,totalValueSB)
                if totalValueSB: staterr[((proc,),(era,),(analysis,analysis+'PP',),(recoSB,))] = 1.+err2SB**0.5/totalValueSB
                results[recoSB]['pp'] = totalValue
                results[recoSB]['ppError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,recoSB) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'PP',),(recoSB,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,recosb)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-totalValueSB)+abs(allShifts[unc+'Down'][mode][mass][recoSB][proc]-totalValueSB))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-allShifts[unc+'Down'][mode][mass][recoSB][proc]))/2.
                        if err and totalValueSB:
                            uncerr[unc][((proc,),(era,),(analysis,analysis+'PP',),(recoSB,))] = min([1+err/totalValueSB,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,reco)] = err/totalValueSB
            # PPR
            for proc in signalsPPR:
                totalValue = 0.
                err2 = 0.
                totalValueSB = 0.
                err2SB = 0.
                for gen in genRecoMap:
                    if len(gen)!=4: continue # 3 for AP, 4 for PP
                    if reco not in genRecoMap[gen]: continue
                    value,err = getCount(counters,proc,'new/allMassWindow/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    scale = scales[mode].scale_Hpp4l(gen[:2],gen[2:])
                    totalValue += scale*value
                    err2 += (scale*err)**2
                    # sb values
                    value,err = getCount(counters,proc,'new/allSideband/{0}/{1}/{2}/gen_{3}'.format(mass,hpphm,reco,gen))
                    totalValueSB += scale*value
                    err2SB += (scale*err)**2
                limits.setExpected(proc,era,analysis+'PPR',reco,totalValue)
                if totalValue: staterr[((proc,),(era,),(analysis+'PPR',),(reco,))] = 1.+err2**0.5/totalValue
                results[reco]['ppR'] = totalValue
                results[reco]['ppRError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,reco) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis+'PPR',),(reco,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,reco)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-totalValue)+abs(allShifts[unc+'Down'][mode][mass][reco][proc]-totalValue))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][reco][proc]-allShifts[unc+'Down'][mode][mass][reco][proc]))/2.
                        if err and totalValue:
                            uncerr[unc][((proc,),(era,),(analysis+'PPR',),(reco,))] = min([1+err/totalValue,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,reco)] = err/totalValue
                # sideband components
                limits.setExpected(proc,era,analysis+'PPR',recoSB,totalValueSB)
                if totalValueSB: staterr[((proc,),(era,),(analysis+'PPR',),(recoSB,))] = 1.+err2SB**0.5/totalValueSB
                results[recoSB]['ppR'] = totalValue
                results[recoSB]['ppRError'] = err2**0.5
                # get the shifts effect on signal
                for unc in shiftTypes:
                    if readUncerr:
                        if '{0}_{1}'.format(proc,recoSB) in uncerr_store[unc]:
                            uncerr[unc][((proc,),(era,),(analysis+'PPR',),(recoSB,))] = 1+uncerr_store[unc]['{0}_{1}'.format(proc,recosb)]
                    else:
                        #err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-totalValueSB)+abs(allShifts[unc+'Down'][mode][mass][recoSB][proc]-totalValueSB))/2.
                        err = (abs(allShifts[unc+'Up'][mode][mass][recoSB][proc]-allShifts[unc+'Down'][mode][mass][recoSB][proc]))/2.
                        if err and totalValueSB:
                            uncerr[unc][((proc,),(era,),(analysis+'PPR',),(recoSB,))] = min([1+err/totalValueSB,2])
                            uncerr_store[unc]['{0}_{1}'.format(proc,reco)] = err/totalValueSB
            # observed
            obs = getCount(counters,'data','new/allMassWindow/{0}/{1}/{2}'.format(mass,hpphm,reco))
            limits.setObserved(era,analysis,reco,obs[0])
            limits.setObserved(era,analysis+'AP',reco,obs[0])
            limits.setObserved(era,analysis+'PP',reco,obs[0])
            limits.setObserved(era,analysis+'PPR',reco,obs[0])
            results[reco]['observed'] = obs[0]
            # sideband
            obs = getCount(counters,'data','new/allSideband/{0}/{1}/{2}'.format(mass,hpphm,reco))
            limits.setObserved(era,analysis,recoSB,obs[0])
            limits.setObserved(era,analysis+'AP',recoSB,obs[0])
            limits.setObserved(era,analysis+'PP',recoSB,obs[0])
            limits.setObserved(era,analysis+'PPR',recoSB,obs[0])
            results[recoSB]['observed'] = obs[0]
            dumpResults(results,'Hpp3l','{0}/{1}'.format(mode,mass))

        # dump the shift uncertainties
        if not readUncerr: dumpResults(uncerr_store,analysis,'shiftUncertainties_{0}_{1}'.format(mode,mass))

        # systematics
        addUncertainties(limits,staterr,uncerr,recoChans,signalsAP+signalsPP+signalsPPR,backgrounds,3)

        # print the datacard
        directory = 'datacards/{0}/{1}'.format('Hpp3l',mode)
        python_mkdir(directory)
        limits.printCard('{0}/{1}'.format(directory,mass),analyses=['Hpp3l'],processes=signalsAP+signalsPP+backgrounds,blind=blind)
        limits.printCard('{0}/{1}AP'.format(directory,mass),analyses=['Hpp3lAP'],processes=signalsAP+backgrounds,blind=blind)
        limits.printCard('{0}/{1}PP'.format(directory,mass),analyses=['Hpp3lPP'],processes=signalsPP+backgrounds,blind=blind)
        limits.printCard('{0}/{1}PPR'.format(directory,mass),analyses=['Hpp3lPPR'],processes=signalsPPR+backgrounds,blind=blind)
