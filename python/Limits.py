import os
import sys
import logging
import numbers

import ROOT

from DevTools.Limits.Models import Model, ModelSpline

class Limits(object):
    '''
    Limits

    A class to encapsulate an analysis selection and produce 
    a datacard that can be read by the Higgs Combine tool.
    '''

    def __init__(self,name='w'):
        self.eras = []        # 7, 8, 13 TeV
        self.analyses = []    # analysis name
        self.channels = []    # analysis channel
        self.observed = {}    # there is one observable per era/analysis/channel combination
        self.processes = {}   # background and signal processes
        self.groups = {}      # groups of systematics
        self.signals = []
        self.backgrounds = []
        self.models = {}      # models to add
        self.expected = {}    # expected yield, one per process/era/analysis/chanel combination
        self.systematics = {} # systematic uncertainties
        self.name = name
        self.workspace = ROOT.RooWorkspace(self.name)

    def __wsimport(self, *args) :
        # getattr since import is special in python
        # NB RooWorkspace clones object
        if len(args) < 2 :
            # Useless RooCmdArg: https://sft.its.cern.ch/jira/browse/ROOT-6785
            args += (ROOT.RooCmdArg(),)
        return getattr(self.workspace, 'import')(*args)

    def __unwrap(self,hist):
        '''Convert 2D histogram to 1D'''
        nbins = hist.GetNbinsX()*hist.GetNbinsY()
        result = ROOT.TH1F(hist.GetName(),hist.GetTitle(),nbins,0,nbins)
        result.SetBinContent(0,hist.GetBinContent(0))
        result.SetBinError(0,hist.GetBinError(0))
        result.SetBinContent(nbins+1,hist.GetBinContent(nbins+1))
        result.SetBinError(nbins+1,hist.GetBinError(nbins+1))
        for nx in range(hist.GetNbinsX()):
            for ny in range(hist.GetNbinsY()):
                b = hist.GetBin(nx+1,ny+1)
                result.SetBinContent(b,hist.GetBinContent(b))
                result.SetBinError(b,hist.GetBinError(b))
        return result

    def addMH(self,mhMin,mhMax):
        self.workspace.factory('MH[{0}, {1}]'.format(mhMin,mhMax))

    def addX(self, xMin, xMax, unit='', label=''):
        self.workspace.factory('x[{0}, {1}]'.format(xMin,xMax))
        if unit: self.workspace.var('x').setUnit(unit)
        if label: self.workspace.var('x').setPlotLabel(label)
        if label: self.workspace.var('x').SetTitle(label)

    def addModel(self,model,label):
        if label in self.models: return
        self.models[label] = model
        self.models[label].build(self.workspace,label)

    def __check(self,test,stored,name='Object'):
        goodToAdd = True
        if 'all' in test: return goodToAdd
        for t in test:
            if t not in stored:
                logging.warning('{0} {1} not recognized.'.format(name,t))
                goodToAdd = False
        return goodToAdd

    def __checkEras(self,eras):
        return self.__check(eras,self.eras,name='Era')

    def __checkAnalyses(self,analyses):
        return self.__check(analyses,self.analyses,name='Analysis')

    def __checkChannels(self,channels):
        return self.__check(channels,self.channels,name='Channel')

    def __checkProcesses(self,processes):
        return self.__check(processes,self.processes,name='Process')

    def addEra(self,era):
        '''Add era to analysis (i.e., 7TeV, 8TeV, 13TeV).'''
        if era in self.eras:
            logging.warning('Era {0} already added.'.format(era))
        else:
            self.eras += [era]

    def addAnalysis(self,analysis):
        '''Add analysis.'''
        if analysis in self.analyses:
            logging.warning('Analysis {0} already added.'.format(analysis))
        else:
            self.analyses += [analysis]

    def addChannel(self,channel):
        '''Add channel to analysis.'''
        if channel in self.channels:
            logging.warning('Channel {0} already added.'.format(channel))
        else:
            self.channels += [channel]

    def addProcess(self,proc,signal=False):
        '''
        Add a process to the datacard.
            parameters:
                signal   - bool            - Declares process to be a signal, default False
        '''
        if proc in self.processes:
            logging.warning('Process {0} already added.'.format(proc))
        else:
            self.processes[proc] = {
                'signal'  : signal,
            }
            if signal:
                self.signals += [proc]
            else:
                self.backgrounds += [proc]

    def addSystematic(self,systname,mode,systematics={}):
        '''
        Add a systematic uncertainty. The name can include the following string
        formatting replacements to set uncorrelated:
            {process}  : uncorrelate process
            {era}      : uncorrelate era
            {analysis} : uncorrelate analysis
            {channel}  : uncorrelate channel
        The supported modes are:
            'lnN'  : log normal uncertainty shape
            'gmN X': gamma function uncertainty shape
            'shape': for shape based uncertainties
            'param': for parameter on unbinned uncertainty
            'flatParam': for freely floating parameter
        The values are set with the 'systematics' arguments. They are dictionaries with the form:
            systematics = {
               (processes,eras,analyses,channels) : value,
            }
        where the key is a tuple of process, era, analysis, and channel the systematic covers, each
        of which is another tuple of the components this sytematic covers.
        'value' is either a number for a rate systematic or a TH1 histogram for a shape uncertainty.
        For asymmetric uncertainties, a tuple should be passed with the first the shift up
        and the second the shift down.
        '''
        if systname in self.systematics:
            logging.warning('Systematic {0} already added.'.format(systname))
        else:
            goodToAdd = True
            for syst in systematics:
                processes,eras,analyses,channels = syst
                goodToAdd = goodToAdd and self.__checkProcesses(processes)
                goodToAdd = goodToAdd and self.__checkEras(eras)
                goodToAdd = goodToAdd and self.__checkAnalyses(analyses)
                goodToAdd = goodToAdd and self.__checkChannels(channels)
            if goodToAdd:
                self.systematics[systname] = {
                    'mode'  : mode,
                    'values': systematics,
                }

    def addGroup(self,groupname,*systnames):
        '''Add a group name for a list of systematics'''
        self.groups[groupname] = systnames

    def getSystematic(self,systname,process,era,analysis,channel):
        '''Return the systematic value for a given systematic/process/era/analysis/channel combination.'''
        # make sure it exists:
        result = 1.
        for syst in self.systematics:
            fullSystName = syst.format(process=process,era=era,analysis=analysis,channel=channel)
            if fullSystName != systname: continue
            # check if there is a systematic value for this combination
            for syst_vals in self.systematics[syst]['values']:
                s_processes, s_eras, s_analyses, s_channels = syst_vals
                if process not in s_processes and 'all' not in s_processes: continue
                if era not in s_eras and 'all' not in s_eras: continue
                if analysis not in s_analyses and 'all' not in s_analyses: continue
                if channel not in s_channels and 'all' not in s_channels: continue
                # return the value
                result = self.systematics[syst]['values'][syst_vals]
        if isinstance(result,ROOT.TH2):
            result = self.__unwrap(result)
        if isinstance(result,tuple) or isinstance(result,list):
            if len(result)==2 and isinstance(result[0],ROOT.TH2):
                result = (self.__unwrap(result[0]),self.__unwrap(result[1]),)
        return result

    def __getSystematicRows(self,syst,processes,era,analysis,channel):
        '''
        Return a dictionary of the systematic values of the form:
            systs = {
                systname : {
                    'mode' : 'lnN',
                    'systs': {
                        (era,analysis,channel,'proc_1'): 1.,
                        ...
                    },
                },
            }
        '''
        systs = {}
        for process in processes:
            systname = syst.format(process=process,era=era,analysis=analysis,channel=channel)
            if systname not in systs:
                systs[systname] = {
                    'mode' : self.systematics[syst]['mode'], 
                    'systs': {},
                }
            systs[systname]['systs'][(era,analysis,channel,process)] = self.getSystematic(systname,process,era,analysis,channel)
        return systs

    def __combineSystematics(self,*systs):
        '''Combine systematics of the form output by __getSystematicRows.'''
        combinedSyst = {}
        for syst in systs:
            for systname in syst:
                if systname not in combinedSyst:
                    combinedSyst[systname] = {
                        'mode' : syst[systname]['mode'],
                        'systs': {},
                    }
                combinedSyst[systname]['systs'].update(syst[systname]['systs'])
        return combinedSyst

    def setObserved(self,era,analysis,channel,value):
        '''Set the observed value for a given era,analysis,channel.'''
        goodToAdd = True
        goodToAdd = goodToAdd and self.__checkEras([era])
        goodToAdd = goodToAdd and self.__checkAnalyses([analysis])
        goodToAdd = goodToAdd and self.__checkChannels([channel])
        if goodToAdd:
            self.observed[(era,analysis,channel)] = value

    def getObserved(self,era,analysis,channel,blind=True,addSignal=False):
        '''Get the observed value. If blinded returns the sum of the expected background.'''
        result = 0.
        if blind:
            bgs = self.backgrounds
            if addSignal: bgs += self.signals
            exp = [self.getExpected(process,era,analysis,channel) for process in bgs]
            if len(exp) and isinstance(exp[0],ROOT.TH1):
                hists = ROOT.TList()
                for e in exp:
                    hists.Add(e)
                if hists.IsEmpty():
                    result = 0.
                else:
                    hist = hists[0].Clone('h_exp')
                    hist.Reset()
                    hist.Merge(hists)
                    for b in range(hist.GetNbinsX()+1):
                        val = int(hist.GetBinContent(b))
                        err = val**0.5
                        hist.SetBinContent(b,val)
                        hist.SetBinError(b,err)
                    result = hist
            else:
                result = sum(exp)
        else:
            key = (era,analysis,channel)
            result = self.observed[key] if key in self.observed else 0.
        if isinstance(result,ROOT.TH2):
            result = self.__unwrap(result)
        return result


    def setExpected(self,process,era,analysis,channel,value):
        '''Set the expected value for a given process,era,analysis,channel.'''
        goodToAdd = True
        goodToAdd = goodToAdd and self.__checkProcesses([process])
        goodToAdd = goodToAdd and self.__checkEras([era])
        goodToAdd = goodToAdd and self.__checkAnalyses([analysis])
        goodToAdd = goodToAdd and self.__checkChannels([channel])
        if goodToAdd:
            self.expected[(process,era,analysis,channel)] = value

    def getExpected(self,process,era,analysis,channel):
        '''Get the expected value.'''
        key = (process,era,analysis,channel)
        val = self.expected[key] if key in self.expected else 0.
        if isinstance(val,ROOT.TH2):
            val = self.__unwrap(val)
        return val if val else 1.0e-10

    def printCard(self,filename,eras=['all'],analyses=['all'],channels=['all'],processes=['all'],blind=True,addSignal=False,saveWorkspace=False,suffix=''):
        '''
        Print a datacard to file.
        Select the eras, analyses, channels you want to include.
        Each will correspond to one bin in the datacard.
        '''

        shapes = self._printMultipleCards(filename,eras,analyses,channels,processes,blind,addSignal,saveWorkspace,suffix)
        
        # shape file
        if shapes:
            outname = filename+'.root'
            if saveWorkspace:
                self.workspace.Print()
                self.workspace.SaveAs(outname)
            else:
                outfile = ROOT.TFile.Open(outname,'RECREATE')
                for shape in shapes:
                    shape.Write('',ROOT.TObject.kOverwrite)
                outfile.Write()
                outfile.Close()

    def _printMultipleCards(self,filename,eras,analyses,channels,processes,blind,addSignal,saveWorkspace,suffix):
        shapes = []
        if isinstance(eras,dict):
            for k,v in eras.iteritems():
                shapes += self._printMultipleCards(filename,v,analyses,channels,processes,blind,addSignal,saveWorkspace,'{0}_{1}'.format(suffix,k))
        elif isinstance(analyses,dict):
            for k,v in analyses.iteritems():
                shapes += self._printMultipleCards(filename,eras,v,channels,processes,blind,addSignal,saveWorkspace,'{0}_{1}'.format(suffix,k))
        elif isinstance(channels,dict):
            for k,v in channels.iteritems():
                shapes += self._printMultipleCards(filename,eras,analyses,v,processes,blind,addSignal,saveWorkspace,'{0}_{1}'.format(suffix,k))
        elif isinstance(processes,dict):
            for k,v in processes.iteritems():
                shapes += self._printMultipleCards(filename,eras,analyses,channels,v,blind,addSignal,saveWorkspace,'{0}_{1}'.format(suffix,k))
        else:
            shapes += self._printSingleCard(filename,eras,analyses,channels,processes,blind,addSignal,saveWorkspace,suffix)
            
        return shapes




    def _printSingleCard(self,filename,eras,analyses,channels,processes,blind,addSignal,saveWorkspace,suffix):
        goodToPrint = True
        goodToPrint = goodToPrint and self.__checkEras(eras)
        goodToPrint = goodToPrint and self.__checkAnalyses(analyses)
        goodToPrint = goodToPrint and self.__checkChannels(channels)
        if not goodToPrint: return

        if eras==['all']: eras = self.eras
        if analyses==['all']: analyses = self.analyses
        if channels==['all']: channels = self.channels
        if processes==['all']: processes = self.processes.keys()
        signals = [x for x in self.signals if x in processes]
        backgrounds = [x for x in self.backgrounds if x in processes]
        shapes = []

        # setup bins
        bins = ['bin']
        observations = ['observation']
        binName = '{era}_{analysis}_{channel}'
        for era in eras:
            for analysis in analyses:
                for channel in channels:
                    blabel = binName.format(era=era,analysis=analysis,channel=channel)
                    bins += [blabel]
                    obs = self.getObserved(era,analysis,channel,blind=blind,addSignal=addSignal)
                    label = 'data_obs_{0}'.format(blabel)
                    if isinstance(obs,ROOT.TH1):
                        logging.debug('{0}: {1}'.format(label,obs.Integral()))
                        obs.SetName(label)
                        obs.SetTitle(label)
                        shapes += [obs]
                        if saveWorkspace:
                            datahist = ROOT.RooDataHist(label, label, ROOT.RooArgList(self.workspace.var("x")), obs)
                            self.__wsimport(datahist)
                            obs = -1
                        else:
                            obs = obs.Integral()
                    else:
                        logging.debug('{0}: {1}'.format(label,obs))
                    # TODO: unbinned data handling
                    observations += ['{0}'.format(obs)]
        imax = len(bins)-1

        # setup processes
        jmax = len(processes)-1

        totalColumns = len(eras)*len(analyses)*len(channels)*len(processes)
        processesOrdered = signals + backgrounds
        binsForRates = ['bin','']+['']*totalColumns
        processNames = ['process','']+['']*totalColumns
        processNumbers = ['process','']+['']*totalColumns
        rates = ['rate','']+['']*totalColumns
        norms = []
        colpos = 1
        for era in eras:
            for analysis in analyses:
                for channel in channels:
                    for process in processesOrdered:
                        colpos += 1
                        binsForRates[colpos] = binName.format(era=era,analysis=analysis,channel=channel)
                        processNames[colpos] = process
                        processNumbers[colpos] = '{0:<10}'.format(processesOrdered.index(process)-len(signals)+1)
                        exp = self.getExpected(process,era,analysis,channel)
                        label = '{0}_{1}'.format(processNames[colpos],binsForRates[colpos])
                        if isinstance(exp,ROOT.TH1):
                            logging.debug('{0}: {1}'.format(label,exp.Integral()))
                            exp.SetName(label)
                            exp.SetTitle(label)
                            shapes += [exp]
                            if saveWorkspace:
                                datahist = ROOT.RooDataHist(label, label, ROOT.RooArgList(self.workspace.var("x")), exp)
                                self.__wsimport(datahist)
                                exp = -1
                            else:
                                exp = exp.Integral()
                        elif isinstance(exp,Model):
                            exp.build(self.workspace,label)
                            if isinstance(exp,ModelSpline):
                                exp = exp.getIntegral(self.workspace)
                            else:
                                exp = exp.getIntegral()
                            norms += [['{}_norm'.format(label),'rateParam',binsForRates[colpos],process,'1']]
                        elif isinstance(exp,numbers.Number):
                            logging.debug('{0}: {1}'.format(label,exp))
                        else:
                            logging.error('Failed to understand: {} {} {} {}'.format(era,analysis,channel,process))
                            print exp
                            raise
                        # TODO: unbinned handling
                        rates[colpos] = '{0:<10.4g}'.format(exp)

        # setup nuissances
        logging.debug('Systs available: {0}'.format([str(x) for x in sorted(self.systematics.keys())]))
        systs = {}
        keys = []
        for era in eras:
            for analysis in analyses:
                for channel in channels:
                    key = (era,analysis,channel)
                    keys += [key]
                    systs[key] = {}
                    for syst in self.systematics:
                        systs[key].update(self.__getSystematicRows(syst,processes,era,analysis,channel))


        combinedSysts = self.__combineSystematics(*[systs[key] for key in systs])
        logging.debug('Systs to add: {0}'.format([str(x) for x in sorted(combinedSysts.keys())]))
        systRows = []
        paramRows = []
        for syst in sorted(combinedSysts.keys()):
            if combinedSysts[syst]['mode'] in ['param','flatParam']:
                # TODO: implement param uncertainties
                pass
            else:
                thisRow = [syst,combinedSysts[syst]['mode']]
                for era in eras:
                    for analysis in analyses:
                        for channel in channels:
                            for process in processesOrdered:
                                key = (era,analysis,channel,process)
                                s = '-'
                                if key in combinedSysts[syst]['systs']:
                                    s = combinedSysts[syst]['systs'][key]
                                    if s==1:
                                        s = '-'
                                    elif isinstance(s,ROOT.TH1):
                                        label = '{0}_{1}_{2}'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                        s.SetName(label)
                                        s.SetTitle(label)
                                        shapes += [s]
                                        if saveWorkspace:
                                            datahist = ROOT.RooDataHist(label, label, ROOT.RooArgList(self.workspace.var("x")), s)
                                            self.__wsimport(datahist)
                                        s = '1'
                                    elif isinstance(s,Model):
                                        label = '{0}_{1}_{2}'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                        s.build(self.workspace,label)
                                        s = '1'
                                    elif (isinstance(s,tuple) or isinstance(s,list)) and len(s)==2:
                                        if isinstance(s[0],ROOT.TH1):
                                            label_up = '{0}_{1}_{2}Up'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                            label_down = '{0}_{1}_{2}Down'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                            s[0].SetName(label_up)
                                            s[0].SetTitle(label_up)
                                            s[1].SetName(label_down)
                                            s[1].SetTitle(label_down)
                                            shapes += s
                                            if saveWorkspace:
                                                datahist_up = ROOT.RooDataHist(label_up, label_up, ROOT.RooArgList(self.workspace.var("x")), s[0])
                                                datahist_down = ROOT.RooDataHist(label_down, label_down, ROOT.RooArgList(self.workspace.var("x")), s[1])
                                                self.__wsimport(datahist_up)
                                                self.__wsimport(datahist_down)
                                            s = '1'
                                        elif isinstance(s[0],Model):
                                            label_up = '{0}_{1}_{2}Up'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                            label_down = '{0}_{1}_{2}Down'.format(process,binName.format(era=era,analysis=analysis,channel=channel),syst)
                                            s[0].build(self.workspace,label_up)
                                            s[1].build(self.workspace,label_down)
                                            s = '1'
                                        elif isinstance(s[0],numbers.Number):
                                            s = '{0:>4.4g}/{1:<4.4g}'.format(*s)
                                        else:
                                            logging.error('Do not know how to handle {0}'.format(s))
                                            raise
                                    elif isinstance(s,numbers.Number):
                                        s = '{0:<10.4g}'.format(s)
                                thisRow += [s]
                systRows += [thisRow]

        kmax = len(systRows)

        # now write to file
        logging.info('Writing {0}{1}.txt'.format(filename,suffix))
        with open(filename+suffix+'.txt','w') as f:
            lineWidth = 80
            firstWidth = 40
            restWidth = 30
            def getline(row):
                try:
                    return '{0} {1}\n'.format(row[0][:firstWidth]+' '*max(0,firstWidth-len(row[0])), ''.join([r[:restWidth]+' '*max(0,restWidth-len(r)) for r in row[1:]]))
                except:
                    print row
                    e = sys.exc_info()[0]
                    print e
                    raise

            def getparamline(row):
                try:
                    return ' '.join([str(x) for x in row])+'\n'
                except:
                    print row
                    e = sys.exc_info()[0]
                    print e
                    raise

            # header
            f.write('imax {0} number of bins\n'.format(imax))
            #f.write('jmax {0} number of processes\n'.format(jmax))
            f.write('jmax * number of processes\n')
            f.write('kmax * number of nuissances\n')
            f.write('-'*lineWidth+'\n')

            # shape information
            if shapes:
                for b in bins[1:]:
                    procString = '$PROCESS_{0}'.format(b)
                    if saveWorkspace: procString = '{0}:{1}'.format(self.name,procString)
                    f.write('shapes * {0} {1}.root {2} {2}_$SYSTEMATIC\n'.format(b,filename,procString))
            else:
                f.write('shapes * * FAKE\n')
            f.write('-'*lineWidth+'\n')
            
            # observation
            f.write(getline(bins))
            f.write(getline(observations))
            f.write('-'*lineWidth+'\n')

            # process definition
            logging.debug('Bins: {0}'.format([str(x) for x in binsForRates]))
            f.write(getline(binsForRates))
            f.write(getline(processNames))
            f.write(getline(processNumbers))
            logging.debug('Rates: {0}'.format([str(x) for x in rates]))
            f.write(getline(rates))
            f.write('-'*lineWidth+'\n')

            # nuissances
            for systRow in systRows:
                logging.debug('Systematic row: {0}'.format([str(x) for x in systRow]))
                f.write(getline(systRow))
            f.write('-'*lineWidth+'\n')

            # rateParams
            for norm in norms:
                logging.debug('Rate param: {0}'.format([str(x) for x in norm]))
                f.write(getparamline(norm))

            # nuissance categories
            for group in self.groups:
                f.write('{0} group = {1}'.format(group,' '.join(self.groups[group])))

        return shapes
