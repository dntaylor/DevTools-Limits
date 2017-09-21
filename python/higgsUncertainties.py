import math

def addUncertainties(limits,staterr,uncerr,recoChans,signals,backgrounds,nl):
    '''Add common uncertainties for the H++ analysis'''
    # stat errs
    limits.addSystematic('stat_{process}_{channel}','lnN',systematics=staterr)

    # shifts
    for unc,errs in uncerr.iteritems():
        limits.addSystematic('{0}_unc'.format(unc),'lnN',systematics=errs)

    systproc = tuple([proc for proc in signals + backgrounds if 'datadriven' not in proc])
    sigproc = tuple([proc for proc in signals])
    approc = tuple([proc for proc in signals if 'HppHmm' not in proc])
    ppproc = tuple([proc for proc in signals if 'HppHmm' in proc])
    ddproc = tuple([proc for proc in backgrounds if 'datadriven' in proc])

    ############
    ### Lumi ###
    ############
    # lumi 2.3% for 2015 and 2.5% for 2016
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
    lumisyst = {
        (systproc,('13TeV80X',),('all',),('all',)): 1.025,
        #(systproc,('13TeV76X',),('all',),('all',)): 1.023,
    }
    limits.addSystematic('lumi','lnN',systematics=lumisyst)

    ################
    ### Electron ###
    ################
    # electron id 2%/leg
    # https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2
    elecsyst = {}
    for c in range(nl):
        systChans = tuple([chan for chan in recoChans if chan.count('e')==c+1]+[chan+'_SB' for chan in recoChans if chan.count('e')==c+1])
        if not systChans: continue
        elecsyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*0.02**2)
    if elecsyst: limits.addSystematic('elec_id','lnN',systematics=elecsyst)

    # electron charge misid 4%
    # TODO update
    #eleccharge = {}
    #for c in range(nl):
    #    systChans = tuple([chan for chan in recoChans if chan.count('e')==c+1]+[chan+'_SB' for chan in recoChans if chan.count('e')==c+1])
    #    if not systChans: continue
    #    eleccharge[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.04**2))
    #if eleccharge: limits.addSystematic('elec_charge','lnN',systematics=eleccharge)

    ############
    ### Muon ###
    ############
    # muon id 1+1%/leg
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
    muonsyst = {}
    for c in range(nl):
        systChans = tuple([chan for chan in recoChans if chan.count('m')==c+1]+[chan+'_SB' for chan in recoChans if chan.count('m')==c+1])
        if not systChans: continue
        muonsyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.01**2 + 0.01**2))
    if muonsyst: limits.addSystematic('muon_id','lnN',systematics=muonsyst)

    # muon single trigger 0.5%
    muontrigsyst = {}
    systChans = tuple([chan for chan in recoChans if chan.count('m')>=1]+[chan+'_SB' for chan in recoChans if chan.count('m')>=1])
    if systChans: muonsyst[(systproc,('all',),('all',),systChans)] = 1.005
    if muontrigsyst: limits.addSystematic('muon_single_trig','lnN',systematics=muontrigsyst)

    ###########
    ### Tau ###
    ###########
    # taus id 5% (SF 0.95)
    # https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV
    tausyst = {}
    for c in range(nl):
        systChans = tuple([chan for chan in recoChans if chan.count('t')==c+1]+[chan+'_SB' for chan in recoChans if chan.count('t')==c+1])
        if not systChans: continue
        tausyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.05**2))
    if tausyst: limits.addSystematic('tau_id','lnN',systematics=tausyst)

    # taus charge misid 2% (SF 1)
    taucharge = {}
    for c in range(nl):
        systChans = tuple([chan for chan in recoChans if chan.count('t')==c+1]+[chan+'_SB' for chan in recoChans if chan.count('t')==c+1])
        if not systChans: continue
        taucharge[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.02**2))
    if taucharge: limits.addSystematic('tau_charge','lnN',systematics=taucharge)

    ##############
    ### Signal ###
    ##############
    # signal unc 15%
    apsyst = {
        (approc, ('all',), ('all',), ('all',)): 1.15,
    }
    ppsyst = {
        (ppproc, ('all',), ('all',), ('all',)): 1.15,
    }
    if approc: limits.addSystematic('sig_unc_AP','lnN',systematics=apsyst)
    if ppproc: limits.addSystematic('sig_unc_PP','lnN',systematics=ppsyst)

    ##################
    ### Datadriven ###
    ##################
    # alpha 10%
    ddsyst = {
        (ddproc, ('all',), ('all',), ('all',)): 1.1,
    }
    if ddproc: limits.addSystematic('alpha_unc','lnN',systematics=ddsyst)

