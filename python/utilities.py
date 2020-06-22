import logging
import os
import sys

import ROOT

def readCount(fileNames,directories,doError=False):
    val = 0.
    err2 = 0.
    for fileName in fileNames:
        if not os.path.isfile(fileName):
            logging.warning('{0} file does not exist.'.format(fileName))
            continue
        tfile = ROOT.TFile.Open(fileName)
        for directory in directories:
            histName = '{0}/count'.format(directory)
            hist = tfile.Get(histName)
            if hist:
                val += hist.GetBinContent(1)
                err2 += hist.GetBinError(1)**2
        tfile.Close()
    return val, err2**0.5 if doError else val

def getDatasetIntegralError(dataset, cutSpec=''):
    
    select = ROOT.RooFormula()
    if cutSpec:
        select = ROOT.RooFormula("select",cutSpec,dataset)
    
    if not cutSpec and not dataset.isWeighted():
        return dataset.numEntries()**0.5
    
    sumw2 = 0
    for i in xrange(dataset.numEntries()):
        dataset.get(i)
        if (cutSpec and select.eval()==0.): continue
        sumw2 += dataset.weight()**2
    
    return sumw2**0.5

def getHistogramIntegralError(hist,binlow=1,binhigh=-1):
    if binhigh<0: binhigh = hist.GetNbinsX()
    integralerr = ROOT.Double(0)
    hist.IntegralAndError(binlow,binhigh,integralerr,"")
    return float(integralerr)
