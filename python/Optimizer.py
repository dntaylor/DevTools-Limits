'''
Tool to optimize selection with several test statistics.
'''
import logging
import itertools

from DevTools.Utilities.utilities import *

class Optimizer(object):
    '''Class to encapsulate the available ways to optimize a selection.'''

    def __init__(self,analysis,**kwargs):
        self.analysis = analysis
        self.statFunc = 'asimovVar'
        self.baseSelection = '1'
        self.mcselection = ''
        self.scalefactor = '1'
        self.mcscalefactor = '1'
        self.datascalefactor = '1'
        self.parameters = {}
        self.signalCounters = []
        self.backgroundCounters = []

    def setTestStatistic(self,statFunc):
        '''Set the test statistic.
           statFunc may be a lambda function or one of the predefined below.
           If it is a function it must take as arguments:
               (0) a string for base selection
               (1+) a selection for this parameter
           Available:
               pois      Poisson (s/sqrt(b))
               poisVar   Poisson with variance in b
               asimov    Asimov
               asimovVar Asimov with variance
        '''
        self.statFunc = statFunc

    def setBaseSelection(self,selection):
        '''Set the base selection to apply before optimization parameters.'''
        self.baseSelection = selection

    def setMCSelection(self,selection):
        '''Set the selection to apply to MC before optimization parameters.'''
        self.mcselection = selection

    def setScalefactor(self,scalefactor):
        '''Set the scalefactor.'''
        self.scalefactor = scalefactor

    def setMCScalefactor(self,scalefactor):
        '''Set the scalefactor.'''
        self.mcscalefactor = scalefactor

    def setDataScalefactor(self,scalefactor):
        '''Set the scalefactor.'''
        self.datascalefactor = scalefactor

    def addParameter(self,varName,var='',logic='',varMin=None,varMax=None,divisions=10,fixed=None):
        '''Add a variable to optimize.
           'logic' is the logic between the variable and the value (>, <, ==)
           Allowed range is set via 'varMin' and 'varMax'.
           The 'divisions' parameter sets the number of values to try in the range.
           The variable can be fixed to a value with the 'fixed' argument.'''
        self.parameters[varName] = {
            'var': var,
            'logic': logic,
            'min': varMin,
            'max': varMax,
            'div': divisions,
            'fix': fixed,
        }

    def addSignalCounters(self,*counters):
        '''Add signal counters (required for available test statistics)'''
        self.signalCounters += counters

    def addBackgroundCounters(self,*counters):
        '''Add background counters (required for available test statistics)'''
        self.backgroundCounters += counters

    def _evaluate(self,func,paramVals):
        return func(self.baseSelection,*paramVals)

    def _getSignalCount(self,selection):
        allCounts = [sumWithError(*c.getCounts('',selection=selection,scalefactor=self.scalefactor,mcscalefactor=self.mcscalefactor,datascalefactor=self.datascalefactor,mccut=self.mcselection).values()) for c in self.signalCounters]
        return sumWithError(*allCounts)

    def _getBackgroundCount(self,selection):
        allCounts = [sumWithError(*c.getCounts('',selection=selection,scalefactor=self.scalefactor,mcscalefactor=self.mcscalefactor,datascalefactor=self.datascalefactor,mccut=self.mcselection).values()) for c in self.backgroundCounters]
        return sumWithError(*allCounts)

    def _getTestStatistic(self,func,*selections):
        selection = ' && '.join(selections)
        sig = self._getSignalCount(selection)
        bg = self._getBackgroundCount(selection)
        return func(sig,bg)

    def _getPrebuilt(self,mode):
        funcMap = {
            'pois'     : poissonSignificance,
            'poisVar'  : poissonSignificanceWithError,
            'asimov'   : asimovSignificance,
            'asimovVar': asimovSignificanceWithError,
        }
        if mode not in funcMap: return None
        return lambda *params: self._getTestStatistic(funcMap[mode],*params)

    def optimize(self):
        '''Run the optimization'''
        func = self._getPrebuilt(self.statFunc)
        if not func: func = self.statFunc
        # setup possible values
        allCuts = {}
        for varName,defs in self.parameters.iteritems():
            if defs['fix']:
                allCuts[varName] = ['{0}{1}{2}'.format(defs['var'],defs['logic'],defs['fix'])]
            elif defs['div'] and defs['min'] and defs['max']:
                allCuts[varName] = ['{0}{1}{2}'.format(defs['var'],defs['logic'],float(x)*(defs['max']-defs['min'])/defs['div']) for x in range(defs['div'])]
            else:
                logging.warning('Invalid parameters for {0}.'.format(varName))
        # iterate through
        testStats = {}
        for params in itertools.product(*[y for x,y in sorted(allCuts.iteritems())]):
            logging.info('Processing: {0}'.format(' '.join(params)))
            testStats[params] = func(self.baseSelection,*params)
            logging.info('          : {0}'.format(testStats[params]))
        return testStats
