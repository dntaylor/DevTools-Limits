import ROOT

from array import array

from DevTools.Plotter.haaUtils import *
import DevTools.Limits.Models as Models


def fitAmm(a,hist,save=False):
    low, high = hist.GetBinLowEdge(1), hist.GetBinUpperEdge(hist.GetNbinsX())
    
    ws = ROOT.RooWorkspace('sig')
    ws.factory('x[{0}, {1}]'.format(low,high))

    model = Models.Voigtian('sig',
        mean  = [a,0,30],
        width = [0.01*a,0,5],
        sigma = [0.01*a,0,5],
    )
    model.build(ws, 'sig')
    results = model.fit(ws, hist, a, save=save)

    return results


def makeSpline(name,masses,xvals,**kwargs):
    save = kwargs.pop('save',False)
    
    ws = ROOT.RooWorkspace('w')
    ws.factory('MH[{0}, {1}]'.format(masses[0],masses[-1]))

    spline = ROOT.RooSpline1D(name, name, ws.var('MH'), len(masses), array('d',masses), array('d',xvals))

    if save:
        numPoints = 100
        mh = ws.var('MH')

        graph = ROOT.TGraph(numPoints)

        orig = mh.getVal()
        for i in range(numPoints):
            p = (masses[-1] - masses[0]) * i/(numPoints-1) + masses[0]
            mh.setVal(p)
            graph.SetPoint(i,p,spline.getVal())
        mh.setVal(orig)

        canvas = ROOT.TCanvas('c','c',800,600)
        graph.Draw()
        graph.SetTitle(name)
        canvas.Print('{0}.png'.format(name))

    return spline


