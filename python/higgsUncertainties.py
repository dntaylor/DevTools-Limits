import math

def addUncertainties(limits,staterr,recoChans,signals,backgrounds):
    '''Add common uncertainties for the H++ analysis'''
    # stat errs
    limits.addSystematic('stat_{process}_{channel}','lnN',systematics=staterr)

    systproc = tuple([proc for proc in signals + backgrounds if proc!='datadriven'])
    sigproc = tuple([proc for proc in signals])
    ddproc = ('datadriven',)

    # lumi 2.7% for 2015 and 6.2% for 2016
    # 2.3% correlated, 5.8% (1.5%) uncorrelated 2016 (2015)
    # https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
    lumisyst = {
        (systproc,('13TeV80X',),('all',),('all',)): 1.023,
        #(systproc,('13TeV76X',),('all',),('all',)): 1.023,
    }
    lumistabilitysyst = {
        (systproc,('13TeV80X',),('all',),('all',)): 1.058,
        #(systproc,('13TeV76X',),('all',),('all',)): 1.015,
    }
    limits.addSystematic('lumi','lnN',systematics=lumisyst)
    limits.addSystematic('lumi_stability_{era}','lnN',systematics=lumistabilitysyst)

    
    # electron id 2%/leg
    elecsyst = {}
    for c in range(3):
        systChans = tuple([chan for chan in recoChans if chan.count('e')==c+1])
        if not systChans: continue
        elecsyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*0.02**2)
    if elecsyst: limits.addSystematic('elec_id','lnN',systematics=elecsyst)

    # muon id 1+0.5%/leg
    muonsyst = {}
    for c in range(3):
        systChans = tuple([chan for chan in recoChans if chan.count('m')==c+1])
        if not systChans: continue
        muonsyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.01**2 + 0.005**2))
    if muonsyst: limits.addSystematic('muon_id','lnN',systematics=muonsyst)

    # muon single trigger 0.5%
    muontrigsyst = {}
    systChans = tuple([chan for chan in recoChans if chan.count('m')>=1])
    if systChans: muonsyst[(systproc,('all',),('all',),systChans)] = 1.005
    if muontrigsyst: limits.addSystematic('muon_single_trig','lnN',systematics=muontrigsyst)

    # taus id 6%
    tausyst = {}
    for c in range(3):
        systChans = tuple([chan for chan in recoChans if chan.count('t')==c+1])
        if not systChans: continue
        tausyst[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.06**2))
    if tausyst: limits.addSystematic('tau_id','lnN',systematics=tausyst)

    # taus charge misid 2.2%
    taucharge = {}
    for c in range(3):
        systChans = tuple([chan for chan in recoChans if chan.count('t')==c+1])
        if not systChans: continue
        taucharge[(systproc,('all',),('all',),systChans)] = 1.+math.sqrt((c+1)*(0.022**2))
    if taucharge: limits.addSystematic('tau_charge','lnN',systematics=taucharge)

    # signal unc 15%
    sigsyst = {
        (sigproc, ('all',), ('all',), ('all',)): 1.15,
    }
    limits.addSystematic('sig_unc','lnN',systematics=sigsyst)

    # alpha 10%
    if 'datadriven' in backgrounds:
        ddsyst = {
            (ddproc, ('all',), ('all',), ('all',)): 1.1,
        }
        limits.addSystematic('alpha_unc_{analysis}_{channel}','lnN',systematics=ddsyst)

