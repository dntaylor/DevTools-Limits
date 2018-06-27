import os
import sys
import logging
import itertools
import numpy as np
import argparse
import math

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch()

from DevTools.Limits.Limits import Limits
from DevTools.Plotter.NtupleWrapper import NtupleWrapper
from DevTools.Utilities.utilities import *
from DevTools.Plotter.haaUtils import *
import DevTools.Limits.Models as Models

logging.basicConfig(level=logging.INFO, stream=sys.stderr, format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

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
    'mm' : [42,4,25],
    'tt' : [60,0,60],
    'h'  : [50,0,1000],
    'hkf': [50,0,1000],
}
rebinning = {
    'mm' : 5, # 10 MeV -> 50 MeV
    'tt' : 1, # 100 MeV -> 100 MeV
    'h'  : 1, # 1 GeV -> 1 GeV
    'hkf': 1, # 1 GeV -> 5 GeV
}

hmasses = [125,300,750]
hmasses = [125]
amasses = [5,7,9,11,13,15,17,19,21]
amasses = [5,11,15,21]
    
signame = 'HToAAH{h}A{a}'
splinename = 'sig{h}'

shiftTypes = ['lep','pu','fake','trig']
shiftTypes = []
shifts = []
for s in shiftTypes:
    shifts += [s+'Up', s+'Down']

#################
### Utilities ###
#################
def getHist(proc,**kwargs):
    scale = kwargs.pop('scale',1)
    shift = kwargs.pop('shift','')
    region = kwargs.pop('region','A')
    do2D = kwargs.pop('do2D',False)
    var = kwargs.pop('var',['mm'])
    wrappers = kwargs.pop('wrappers',{})
    if do2D:
        plot = '{}_{}'.format(*[varHists[v] for v in var])
    else:
        plot = varHists[var[0]]
    plotname = 'region{}/{}'.format(region,plot)
    if do2D:
        hists = [wrappers[s+shift].getHist2D(plotname) for s in sampleMap[proc]]
    else:
        hists = [wrappers[s+shift].getHist(plotname) for s in sampleMap[proc]]
    hist = sumHists(proc+region,*hists)
    hist.Scale(scale)
    return hist

def getDatadrivenHist(**kwargs):
    shift = kwargs.pop('shift','')
    source = kwargs.pop('source','B')
    region = kwargs.pop('region','A')
    do2D = kwargs.pop('do2D',False)
    var = kwargs.pop('var',['mm'])
    wrappers = kwargs.pop('wrappers',{})
    if do2D:
        plot = '{}_{}'.format(*[varHists[v] for v in var])
    else:
        plot = varHists[var[0]]
    plotname = 'region{}_fakeFor{}/{}'.format(source,region,plot)
    if do2D:
        hists = [wrappers[s+shift].getHist2D(plotname) for s in sampleMap['data']]
    else:
        hists = [wrappers[s+shift].getHist(plotname) for s in sampleMap['data']]
    hist = sumHists('data'+region+source,*hists)
    return hist

def getMatrixHist(proc,**kwargs):
    scale = kwargs.pop('scale',1)
    shift = kwargs.pop('shift','')
    region = kwargs.pop('region','A')
    sources = kwargs.pop('sources',['A','C'])
    doPrompt = kwargs.pop('doPrompt',True)
    doFake = kwargs.pop('doFake',False)
    do2D = kwargs.pop('do2D',False)
    var = kwargs.pop('var',['mm'])
    wrappers = kwargs.pop('wrappers',{})
    if do2D:
        plot = '{}_{}'.format(*[varHists[v] for v in var])
    else:
        plot = varHists[var[0]]
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
    hist.Scale(scale)
    return hist

def getMatrixDatadrivenHist(**kwargs):
    shift = kwargs.pop('shift','')
    region = kwargs.pop('region','A')
    sources = kwargs.pop('sources',['A','C'])
    fakeRegion = kwargs.pop('fakeRegion','B')
    fakeSources = kwargs.pop('fakeSources',['B','D'])
    doPrompt = kwargs.pop('doPrompt',True)
    doFake = kwargs.pop('doFake',False)
    do2D = kwargs.pop('do2D',False)
    var = kwargs.pop('var',['mm'])
    wrappers = kwargs.pop('wrappers',{})
    if do2D:
        plot = '{}_{}'.format(*[varHists[v] for v in var])
    else:
        plot = varHists[var[0]]
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

def getSpline(histMap,h,var=['mm'],tag=''):
    # initial fit
    results = {}
    errors = {}
    results[h] = {}
    errors[h] = {}
    for a in amasses:
        ws = ROOT.RooWorkspace('sig')
        binning = varBinning[var[0]]
        ws.factory('x[{0}, {1}]'.format(*binning[1:]))
        ws.var('x').setUnit('GeV')
        ws.var('x').setPlotLabel('m_{#mu#mu}')
        ws.var('x').SetTitle('m_{#mu#mu}')
        model = Models.Voigtian('sig',
            mean  = [a,0,30],
            width = [0.01*a,0,5],
            sigma = [0.01*a,0,5],
        )
        model.build(ws, 'sig')
        hist = histMap[signame.format(h=h,a=a)]
        results[h][a], errors[h][a] = model.fit(ws, hist, '{0}_{1}{2}'.format(h,a,tag), save=True, doErrors=True)

    models = {
        'mean' : Models.Chebychev('mean',  order = 1, p0 = [0,-1,1], p1 = [0.1,-1,1], p2 = [0.03,-1,1]),
        'width': Models.Chebychev('width', order = 1, p0 = [0,-1,1], p1 = [0.1,-1,1], p2 = [0.03,-1,1]),
        'sigma': Models.Chebychev('sigma', order = 1, p0 = [0,-1,1], p1 = [0.1,-1,1], p2 = [0.03,-1,1]),
    }

    for param in ['mean', 'width', 'sigma']:
        ws = ROOT.RooWorkspace(param)
        ws.factory('x[{},{}]'.format(1,30))
        ws.var('x').setUnit('GeV')
        ws.var('x').setPlotLabel('m_{#mu#mu}')
        ws.var('x').SetTitle('m_{#mu#mu}')
        model = models[param]
        model.build(ws, param)
        name = '{}_{}{}'.format(param,h,tag)
        hist = ROOT.TH1D(name, name, len(amasses), 4, 22)
        vals = [results[h][a]['{}_{}_{}{}'.format(param,h,a,tag)] for a in amasses]
        errs = [errors[h][a]['{}_{}_{}{}'.format(param,h,a,tag)] for a in amasses]
        for i,a in enumerate(amasses):
            b = hist.FindBin(a)
            hist.SetBinContent(b,vals[i])
            hist.SetBinError(b,errs[i])
        model.fit(ws, hist, name, save=True)

    # create model
    for a in amasses:
        print h, a, results[h][a]
    model = Models.VoigtianSpline(splinename.format(h=h),
        **{
            'masses' : amasses,
            'means'  : [results[h][a]['mean_{0}_{1}{2}'.format(h,a,tag)] for a in amasses],
            'widths' : [results[h][a]['width_{0}_{1}{2}'.format(h,a,tag)] for a in amasses],
            'sigmas' : [results[h][a]['sigma_{0}_{1}{2}'.format(h,a,tag)] for a in amasses],
        }
    )
    integrals = [histMap[signame.format(h=h,a=a)].Integral() for a in amasses]
    model.setIntegral(amasses,integrals)

    return model
        
def buildModel(ws,**kwargs):
    tag = kwargs.pop('tag','')
    isLimits = isinstance(ws,Limits)


    # upsilon
    upsilon1S = Models.Gaussian('upsilon1S',
        mean  = [9.5,9.3,9.7],
        sigma = [0.1,0,0.3],
    )
    name = 'upsilon1S{}'.format('_'+tag if tag else '')
    if isLimits:
        ws.addModel(upsilon1S,name)
    else:
        upsilon1S.build(ws,name)

    upsilon2S = Models.Gaussian('upsilon2S',
        mean  = [10.0,9.8,10.2],
        sigma = [0.1,0,0.3],
    )
    name = 'upsilon2S{}'.format('_'+tag if tag else '')
    if isLimits:
        ws.addModel(upsilon2S,name)
    else:
        upsilon2S.build(ws,name)

    upsilon3S = Models.Gaussian('upsilon3S',
        mean  = [10.3,10.2,10.4],
        sigma = [0.1,0,0.3],
    )
    name = 'upsilon3S{}'.format('_'+tag if tag else '')
    if isLimits:
        ws.addModel(upsilon3S,name)
    else:
        upsilon3S.build(ws,name)

    upsilon = Models.Sum('upsilon',
        **{
            'upsilon1S{}'.format('_'+tag if tag else '') : [0,1],
            'upsilon2S{}'.format('_'+tag if tag else '') : [0,1],
            'upsilon3S{}'.format('_'+tag if tag else '') : [0,1],
            'recursive' : True,
        }
    )
    name = 'upsilon{}'.format('_'+tag if tag else '')
    if isLimits:
        ws.addModel(upsilon,name)
    else:
        upsilon.build(ws,name)

    # continuum background
    cont = Models.Chebychev('cont',
        order = 2,
        p0 = [-1,-2,0],
        p1 = [0.1,0,0.5],
        p2 = [0.03,-1,1],
    )
    name = 'cont{}'.format('_'+tag if tag else '')
    if isLimits:
        ws.addModel(cont,name)
    else:
        cont.build(ws,name)

    # sum
    bg = Models.Sum('bg',
        **{
            'cont{}'.format('_'+tag if tag else '')    : [0,1],
            'upsilon{}'.format('_'+tag if tag else '') : [0,1],
            'recursive' : True,
        }
    )
    name = 'bg{}'.format('_'+tag if tag else '')
    if isLimits:
        return bg
    else:
        bg.build(ws,name)



###############
### Control ###
###############

def create_datacard(args):
    doMatrix = False
    doParametric = args.parametric
    doUnbinned = args.unbinned
    do2D = len(args.fitVars)==2
    blind = not args.unblind
    addSignal = args.addSignal
    signalParams = {'h': args.higgs, 'a': args.pseudoscalar}
    wsname = 'w'
    var = args.fitVars
    
    if do2D and doParametric:
       logging.error('Parametric 2D fits are not yet supported')
       raise

    if doUnbinned and not doParametric:
        logging.error('Unbinned only supported with parametric option')
        raise
    

    #############
    ### Setup ###
    #############
    sampleMap = getSampleMap()
    
    backgrounds = ['datadriven']
    data = ['data']
    
    signals = [signame.format(h=h,a=a) for h in hmasses for a in amasses]
    signalToAdd = signame.format(**signalParams)
    signalSplines = [splinename.format(h=h) for h in hmasses]

    
    wrappers = {}
    for proc in backgrounds+signals+data:
        if proc=='datadriven': continue
        for sample in sampleMap[proc]:
            wrappers[sample] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X')
            for shift in shifts:
                wrappers[sample+shift] = NtupleWrapper('MuMuTauTau',sample,new=True,version='80X',shift=shift)
    
    ##############################
    ### Create/read histograms ###
    ##############################
    
    histMap = {}
    # The definitons of which regions match to which arguments
    # PP can take a fake rate datadriven estimate from PF, but PF can only take the observed values
    regionArgs = {
        'PP': {'region':'A','fakeRegion':'B','source':'B','sources':['A','C'],'fakeSources':['B','D'],},
        'PF': {'region':'B','sources':['B','D'],},
    }
    for mode in ['PP','PF']:
        histMap[mode] = {}
        for shift in ['']+shifts:
            histMap[mode][shift] = {}
            for proc in backgrounds+signals:
                logging.info('Getting {} {}'.format(proc,shift))
                if proc=='datadriven':
                    # TODO: unbinned, get the RooDataHist from flattenener first
                    if mode=='PP':
                        if doMatrix:
                            histMap[mode][shift][proc] = getMatrixDatadrivenHist(var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                        else:
                            histMap[mode][shift][proc] = getDatadrivenHist(var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                    else:
                        if doMatrix:
                            histMap[mode][shift][proc] = getMatrixHist('data',var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                        else:
                            histMap[mode][shift][proc] = getHist('data',var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                else:
                    if doMatrix:
                        histMap[mode][shift][proc] = getMatrixHist(proc,var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                    else:
                        histMap[mode][shift][proc] = getHist(proc,var=var,wrappers=wrappers,shift=shift,do2D=do2D,**regionArgs[mode])
                if do2D:
                    pass # TODO, figure out how to rebin 2D
                else:
                    histMap[mode][shift][proc].Rebin(rebinning[var[0]])
            if shift: continue
            logging.info('Getting observed')
            if blind:
                samples = backgrounds
                if addSignal: samples = backgrounds + [signalToAdd]
                hists = []
                for proc in samples:
                    hists += [histMap[mode][shift][proc]]
                hist = sumHists('obs',*hists)
                #for b in range(hist.GetNbinsX()+1):
                #    val = int(hist.GetBinContent(b))
                #    if val<0: val = 0
                #    err = val**0.5
                #    hist.SetBinContent(b,val)
                #    #hist.SetBinError(b,err)
                histMap[mode][shift]['data'] = hist
            else:
                hist = getHist('data',var=var,wrappers=wrappers,do2D=do2D,**regionArgs[mode])
                histMap[mode][shift]['data'] = hist
                if do2D:
                    pass
                else:
                    histMap[mode][shift]['data'].Rebin(rebinning[var[0]])
    
    #####################
    ### Create Limits ###
    #####################
    limits = Limits(wsname)
    
    limits.addEra('Run2016')
    limits.addAnalysis('HAA')
    
    era = 'Run2016'
    analysis = 'HAA'
    reco = 'mmmt'
    
    for mode in ['PP','PF']:
        limits.addChannel(mode)
        if doParametric:
            binning = varBinning[var[0]]
            limits.addMH(*binning[1:])
            limits.addX(*binning[1:],unit='GeV',label='m_{#mu#mu}')
            for h in hmasses:
                limits.addProcess(splinename.format(h=h),signal=True)
            for background in backgrounds:
                limits.addProcess(background)
            
            # add models
            for h in hmasses:
                model = getSpline(histMap[mode][''],h,tag=mode)
                limits.setExpected(splinename.format(h=h),era,analysis,mode,model)

            if doUnbinned:
                bg = buildModel(limits,tag=mode)
                limits.setExpected('datadriven', era, analysis, mode, bg)
            else:
                # add histograms for background if not using an unbinned model
                for bg in backgrounds:
                    limits.setExpected(bg,era,analysis,mode,histMap[mode][''][bg])
            
            # get roodatahist
            limits.setObserved(era,analysis,mode,histMap[mode]['']['data'])
        
        else:
        
            for signal in signals:
                limits.addProcess(signal,signal=True)
            for background in backgrounds:
                limits.addProcess(background)
            
            for proc in backgrounds:
                limits.setExpected(proc,era,analysis,mode,histMap[mode][''][proc])
            for proc in signals:
                limits.setExpected(proc,era,analysis,mode,histMap[mode][''][proc])
            
            limits.setObserved(era,analysis,mode,histMap[mode]['']['data'])
        
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
        statMapUp[proc] = getStat(histMap[mode][''][proc],'Up')
        statMapDown[proc] = getStat(histMap[mode][''][proc],'Down')
    statsyst = {}

    for mode in ['PP','PF']:
        # background
        if doUnbinned:
            # TODO: add errors on params
            pass
        else:
            for proc in bgproc:
                statsyst[((proc,),(era,),(analysis,),(mode,))] = (statMapUp[proc],statMapDown[proc])

        # signal
        if doParametric:
            for h in hmasses:
                statsyst[((splinename.format(h=h),),(era,),(analysis,),(mode,))] = (getSpline(statMapUp,h,tag=mode+'StatUp'),getSpline(statMapDown,h,tag=mode+'StatDown'))
        else:
            for proc in sigproc:
                statsyst[((proc,),(era,),(analysis,),(mode,))] = (statMapUp[proc],statMapDown[proc])

    limits.addSystematic('stat_{process}_{channel}','shape',systematics=statsyst)

    ##############
    ### shifts ###
    ##############
    for shift in shiftTypes:
        logging.info('Adding {} systematic'.format(shift))
        shiftsyst = {}
        for mode in ['PP','PF']:

            # background
            if doUnbinned:
                # TODO rateParams on bg model
                pass
            else:
                for proc in bgproc:
                    shiftsyst[((proc,),(era,),(analysis,),(mode,))] = (histMap[mode][shift+'Up'][proc], histMap[mode][shift+'Down'][proc])

            # signal
            if doParametric:
                for h in hmasses:
                    shiftsyst[((splinename.format(h=h),),(era,),(analysis,),(mode,))] = (getSpline(histMap[mode][shift+'Up'],h,tag=mode+shift+'Up'),getSpline(histMap[mode][shift+'Down'],h,tag=mode+shift+'Down'))
            else:
                for proc in sigproc:
                    shiftsyst[((proc,),(era,),(analysis,),(mode,))] = (histMap[mode][shift+'Up'][proc], histMap[mode][shift+'Down'][proc])

        limits.addSystematic(shift,'shape',systematics=shiftsyst)
        
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
        (muproc,(era,),('all',),('all',)): 1+math.sqrt(sum([0.01**2,0.005**2]*2+[0.01**2])), # 2 lead have iso, tau_mu doesnt
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
    parser.add_argument('--unbinned', action='store_true', help='Create unbinned datacards')
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
