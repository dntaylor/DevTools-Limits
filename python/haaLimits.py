import os
import sys
import logging
import itertools
import numpy as np
import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
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

def create_datacard(args):
    baseCut = 'fabs(am1_dxy)<0.2 && fabs(am1_dz)<0.5 && fabs(am2_dxy)<0.2 && fabs(am2_dz)<0.5 && fabs(atm_dxy)<0.2 && fabs(atm_dz)<0.5 && fabs(ath_dz)<0.5'
    regions = {
        'A': '(am1_isolation<0.15 && am2_isolation<0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT>0.5' + ' && ' + baseCut,
        'B': '(am1_isolation<0.15 && am2_isolation<0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT<0.5' + ' && ' + baseCut,
        'C': '(am1_isolation>0.15 || am2_isolation>0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT>0.5' + ' && ' + baseCut,
        'D': '(am1_isolation>0.15 || am2_isolation>0.15) && ath_byVLooseIsolationMVArun2v1DBoldDMwLT<0.5' + ' && ' + baseCut,
    }
    
    doMatrix = True
    doParametric = args.parametric
    do2D = len(args.fitVars)==2
    blind = not args.unblind
    addSignal = args.addSignal
    selection = regions['A']
    scaleVars = ['genWeight','pileupWeight','triggerEfficiency']
    scalefactor = '*'.join(scaleVars)
    signalParams = {'h': args.higgs, 'a': args.pseudoscalar}
    wsname = 'w'
    var = args.fitVars
    
    if do2D and doParametric:
       logging.error('Parametric 2D fits are not yet supported')
       raise
    
    varHists = {
        'mm' : 'ammMass',
        'tt' : 'attMass',
        'h'  : 'hMass',
        'hkf': 'hMassKinFit',
    }
    varNames = {
        'mm' : 'amm_mass',
        'tt' : 'att_mass',
        'h'  : 'h_mass',
        'hkf': 'h_massKinFit',
    }
    varBinning = {
        'mm' : [42,4,25] if do2D else [420,4,25],
        'tt' : [60,0,60] if do2D else [600,0,60],
        'h'  : [50,0,1000] if do2D else [1000,0,1000],
        'hkf': [50,0,1000] if do2D else [1000,0,1000],
    }
    rebinning = {
        'mm' : 5, # 10 MeV -> 50 MeV
        'tt' : 1, # 100 MeV -> 100 MeV
        'h'  : 1, # 1 GeV -> 1 GeV
        'kkf': 1, # 1 GeV -> 5 GeV
    }

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
    splinename = 'sig{h}'
    
    hmasses = [125,300,750]
    #hmasses = [125]
    amasses = [5,7,9,11,13,15,17,19,21]
    #amasses = [5,11,15,21]
    
    signals = [signame.format(h=h,a=a) for h in hmasses for a in amasses]
    signalToAdd = signame.format(**signalParams)
    signalSplines = [splinename.format(h=h) for h in hmasses]

    shifts = ['lepUp','lepDown','puUp','puDown','fakeUp','fakeDown','trigUp','trigDown']
    
    wrappers = {}
    for proc in backgrounds+signals+data:
        if proc=='datadriven': continue
        for sample in sampleMap[proc]:
            wrappers[sample] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X')
            for shift in shifts:
                wrappers[sample+shift] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X',shift=shift)
    
    #################
    ### Utilities ###
    #################
    def getHist(proc,**kwargs):
        scale = kwargs.pop('scale',1)
        shift = kwargs.pop('shift','')
        if do2D:
            plot = '{}_{}'.format(*[varHists[v] for v in var])
        else:
            plot = varHists[var[0]]
        region = 'A'
        plotname = 'region{}/{}'.format(region,plot)
        if do2D:
            hists = [wrappers[s+shift].getHist2D(plotname) for s in sampleMap[proc]]
        else:
            hists = [wrappers[s+shift].getHist(plotname) for s in sampleMap[proc]]
        hist = sumHists(proc+region,*hists)
        #hist.Rebin(10)
        hist.Scale(scale)
        return hist
    
    def getDatadrivenHist(**kwargs):
        shift = kwargs.pop('shift','')
        if do2D:
            plot = '{}_{}'.format(*[varHists[v] for v in var])
        else:
            plot = varHists[var[0]]
        source = 'B'
        region = 'A'
        plotname = 'region{}_fakeFor{}/{}'.format(source,region,plot)
        if do2D:
            hists = [wrappers[s+shift].getHist2D(plotname) for s in sampleMap['data']]
        else:
            hists = [wrappers[s+shift].getHist(plotname) for s in sampleMap['data']]
        hist = sumHists('data'+region+source,*hists)
        #hist.Rebin(10)
        return hist
    
    def getMatrixHist(proc,**kwargs):
        scale = kwargs.pop('scale',1)
        shift = kwargs.pop('shift','')
        if do2D:
            plot = '{}_{}'.format(*[varHists[v] for v in var])
        else:
            plot = varHists[var[0]]
        region = 'A'
        sources = ['A','C']
        fakeRegion = 'B'
        fakeSources = ['B','D']
        doPrompt = True
        doFake = False
        applot = ['matrixP/region{}_for{}/{}'.format(source,region,plot) for source in sources]
        afplot = ['matrixF/region{}_for{}/{}'.format(source,region,plot) for source in sources]
        hists = []
        for s in sampleMap[proc]:
            if do2D:
                if doPrompt: hists += [wrappers[s+shift].getHist2D(plotname) for plotname in applot]
                if doFake: hists += [wrappers[s+shift].getHist2D(plotname) for plotname in afplot]
            else:
                if doPrompt: hists += [wrappers[s+shift].getHist(plotname) for plotname in applot]
                if doFake: hists += [wrappers[s+shift].getHist(plotname) for plotname in afplot]
        hist = sumHists(proc+region+source,*hists)
        #hist.Rebin(10)
        hist.Scale(scale)
        return hist
    
    def getMatrixDatadrivenHist(**kwargs):
        shift = kwargs.pop('shift','')
        if do2D:
            plot = '{}_{}'.format(*[varHists[v] for v in var])
        else:
            plot = varHists[var[0]]
        region = 'A'
        sources = ['A','C']
        fakeRegion = 'B'
        fakeSources = ['B','D']
        doPrompt = True
        doFake = False
        bpplot = ['matrixP/region{}_for{}_fakeFor{}/{}'.format(source,fakeRegion,region,plot) for source in fakeSources]
        bfplot = ['matrixF/region{}_for{}_fakeFor{}/{}'.format(source,fakeRegion,region,plot) for source in fakeSources]
        hists = []
        for s in sampleMap['data']:
            if do2D:
                if doPrompt: hists += [wrappers[s+shift].getHist2D(plotname) for plotname in bpplot]
                if doFake: hists += [wrappers[s+shift].getHist2D(plotname) for plotname in bfplot]
            else:
                if doPrompt: hists += [wrappers[s+shift].getHist(plotname) for plotname in bpplot]
                if doFake: hists += [wrappers[s+shift].getHist(plotname) for plotname in bfplot]
        hist = sumHists('data'+region+source,*hists)
        #hist.Rebin(10)
        return hist
    
    def getBinned(proc,**kwargs):
        sf = kwargs.pop('scalefactor','1' if proc=='data' else scalefactor)
        sel = kwargs.pop('selection',selection)
        shift = kwargs.pop('shift','')
    
        hists = ROOT.TList()
        for sample in sampleMap[proc]:
            if do2D:
                hist = wrappers[sample+shift].getTempHist2D(sample,sel,sf,varNames[var[0]],varNames[var[1]],varBinning[var[0]],varBinning[var[1]])
            else:
                hist = wrappers[sample+shift].getTempHist(sample,sel,sf,varNames[var[0]],varBinning[var[0]])
            hists.Add(hist)
        if hists.IsEmpty():
            hist = 0
        else:
            hist = hists[0].Clone('h_{0}'.format(proc))
            hist.Reset()
            hist.Merge(hists)
        return hist

    def getDatadrivenTemp(**kwargs):
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
            binning = varBinning[var[0]]
            ws.factory('x[{0}, {1}]'.format(*binning[1:]))
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
        model.setIntegral(amasses,integrals)
    
        return model
    
    ##############################
    ### Create/read histograms ###
    ##############################
    
    histMap = {}
    for proc in backgrounds+signals:
        logging.info('Getting {0}'.format(proc))
        if proc=='datadriven':
            #histMap[proc] = getDatadriven()
            if doMatrix:
                histMap[proc] = getMatrixDatadrivenHist()
            else:
                histMap[proc] = getDatadrivenHist()
        else:
            #histMap[proc] = getBinned(proc)
            scale = 0.05 if proc in signals else 1
            if doMatrix:
                histMap[proc] = getMatrixHist(proc,scale=scale)
            else:
                histMap[proc] = getHist(proc,scale=scale)
        if do2D:
            pass # TODO, figure out how to rebin 2D
        else:
            histMap[proc].Rebin(rebinning[var[0]])
    logging.info('Getting observed')
    if blind:
        samples = backgrounds
        if addSignal: samples = backgrounds + [signalToAdd]
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
        binning = varBinning[var[0]]
        limits.addMH(*binning[1:])
        limits.addX(*binning[1:])
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
        # TODO, uncertainty on parameter used to interpolate between
        pass
        #for h in hmasses:
        #    statsyst[((splinename.format(h=h),),(era,),(analysis,),(reco,))] = (getSpline(statMapUp,h),getSpline(statMapDown,h))
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
    lumiproc = systsplineproc if doParametric else systproc
    lumisyst = {
        (lumiproc,(era,),('all',),('all',)): 1.025,
    }
    limits.addSystematic('lumi','lnN',systematics=lumisyst)

    ############
    ### muon ###
    ############
    # from z: 1 % + 0.5 % + 0.5 % per muon for id + iso + trig (pt>20)
    logging.info('Adding mu id+iso systematic')
    muproc = systsplineproc if doParametric else systproc
    musyst = {
        (muproc,(era,),('all',),('all',)): 1+sqrt(sum([0.01**2,0.005**2]*2+[0.01**2])), # 2 lead have iso, tau_mu doesnt
    }
    limits.addSystematic('muid','lnN',systematics=musyst)

    logging.info('Adding mu trig systematic')
    musyst = {
        (muproc,(era,),('all',),('all',)): 1.005, # 1 triggering muon
    }
    limits.addSystematic('mutrig','lnN',systematics=musyst)

    ###########
    ### tau ###
    ###########
    # 5% on sf 0.99 (VL/L) or 0.97 (M)
    logging.info('Adding mu id+iso systematic')
    tauproc = systsplineproc if doParametric else systproc
    tausyst = {
        (tauproc,(era,),('all',),('all',)): 1.05,
    }
    limits.addSystematic('tauid','lnN',systematics=tausyst)

    
    ##############
    ### Pileup ###
    ##############
    #logging.info('Adding pileup systematic')
    #puMapUp = {}
    #puMapDown = {}
    #for proc in backgrounds+signals:
    #    logging.info('Getting {0} PU up'.format(proc))
    #    if proc=='datadriven':
    #        puMapUp[proc] = getDatadriven(scalefactor='*'.join(['pileupWeightUp' if x=='pileupWeight' else x for x in scaleVars]))
    #    else:
    #        puMapUp[proc] = getBinned(proc,scalefactor='*'.join(['pileupWeightUp' if x=='pileupWeight' else x for x in scaleVars]))
    #    logging.info('Getting {0} PU down'.format(proc))
    #    if proc=='datadriven':
    #        puMapDown[proc] = getDatadriven(scalefactor='*'.join(['pileupWeightDown' if x=='pileupWeight' else x for x in scaleVars]))
    #    else:
    #        puMapDown[proc] = getBinned(proc,scalefactor='*'.join(['pileupWeightDown' if x=='pileupWeight' else x for x in scaleVars]))
    #pusyst = {}
    #for proc in bgproc:
    #    pusyst[((proc,),(era,),(analysis,),(reco,))] = (puMapUp[proc],puMapDown[proc])
    #if doParametric:
    #    # TODO, uncertainty on parameter used to interpolate between
    #    pass
    #    #for h in hmasses:
    #    #    pusyst[((splinename.format(h=h),),(era,),(analysis,),(reco,))] = (getSpline(puMapUp,h),getSpline(puMapDown,h))
    #else:
    #    for proc in sigproc:
    #        pusyst[((proc,),(era,),(analysis,),(reco,))] = (puMapUp[proc],puMapDown[proc])
    #limits.addSystematic('pu','shape',systematics=pusyst)
    
    ######################
    ### Print datacard ###
    ######################
    directory = 'datacards_shape/{0}'.format('MuMuTauTau')
    python_mkdir(directory)
    datacard = '{0}/mmmt_{1}'.format(directory, args.tag) if args.tag else '{}/mmmt'.format(directory)
    processes = {}
    if doParametric:
        for h in hmasses:
            processes[signame.format(h=h,a='X')] = [splinename.format(h=h)] + backgrounds
    else:
        for signal in signals:
            processes[signal] = [signal]+backgrounds
    limits.printCard(datacard,processes=processes,blind=False,saveWorkspace=doParametric)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Create datacard')

    parser.add_argument('fitVars', type=str, nargs='*', default=[])
    parser.add_argument('--unblind', action='store_true', help='Unblind the datacards')
    parser.add_argument('--parametric', action='store_true', help='Create parametric datacards')
    parser.add_argument('--addSignal', action='store_true', help='Insert fake signal')
    parser.add_argument('--higgs', type=int, default=125, choices=[125,300,750])
    parser.add_argument('--pseudoscalar', type=int, default=15, choices=[5,7,9,11,13,15,17,19,21])
    parser.add_argument('--tag', type=str, default='')

    return parser.parse_args(argv)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    create_datacard(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
