import logging

from array import array

import ROOT



class Model(object):

    def __init__(self,name,**kwargs):
        self.name = name
        self.x = kwargs.pop('x','x')
        self.y = kwargs.pop('y','y')
        self.z = kwargs.pop('z','z')
        self.kwargs = kwargs

    def update(self,**kwargs):
        '''Update the floating parameters'''
        self.kwargs.update(kwargs)

    def build(self,ws,label):
        '''Dummy method to add model to workspace'''
        pass

    def fit(self,ws,hist,name,save=False):
        '''Fit the model to a histogram and return the fit values'''
        if isinstance(hist,ROOT.TH1):
            dhname = 'dh_{0}'.format(name)
            hist = ROOT.RooDataHist(dhname, dhname, ROOT.RooArgList(ws.var(self.x)), hist)
        self.build(ws,name)
        model = ws.pdf(name)
        fr = model.fitTo(hist,ROOT.RooFit.Save())
        pars = fr.floatParsFinal()
        vals = {}
        for p in range(pars.getSize()):
            vals[pars.at(p).GetName()] = pars.at(p).getValV()

        if save:
            savename = '{0}_{1}'.format(self.name,name)
            x = ws.var(self.x)
            xFrame = x.frame()
            hist.plotOn(xFrame)
            model.plotOn(xFrame)
            model.paramOn(xFrame)
            canvas = ROOT.TCanvas(savename,savename,800,800)
            xFrame.Draw()
            canvas.Print('{0}.png'.format(savename))

        return vals

    def setIntegral(self,integral):
        self.integral = integral

    def getIntegral(self):
        if hasattr(self,'integral'): return self.integral
        return 1

class ModelSpline(Model):

    def __init__(self,name,**kwargs):
        self.MH = kwargs.pop('MH','MH')
        super(ModelSpline,self).__init__(name,**kwargs)

    def setIntegral(self,masses,integrals):
        self.masses = masses
        self.integral = integrals

    def getIntegral(self,ws):
        if not hasattr(self,'integral'): return 1
        integralName = '{0}_norm'.format(self.name) # TODO, better name
        integralSpline  = ROOT.RooSpline1D(integralName,  integralName,  ws.var(self.MH), len(self.masses), array('d',self.masses), array('d',self.integrals))
        # import to workspace
        getattr(ws, "import")(integralSpline, ROOT.RooFit.RecycleConflictNodes())
        # return name
        return integralName

class SplineParam(object):

    def __init__(self,name,**kwargs):
        self.name = name
        self.kwargs = kwargs

    def build(self,ws,label):
        paramName = '{0}'.format(label) 
        ws.factory('{0}[0,-10,10]'.format(paramName))

class Spline(object):

    def __init__(self,name,**kwargs):
        self.name = name
        self.mh = kwargs.pop('MH','MH')
        self.kwargs = kwargs

    def build(self,ws,label):
        masses = self.kwargs.get('masses', [])
        values = self.kwargs.get('values', [])
        shifts = self.kwargs.get('shifts', {})
        splineName = label
        if shifts:
            centralName = '{0}_central'.format(label)
            splineCentral = ROOT.RooSpline1D(centralName,  centralName,  ws.var(self.mh), len(masses), array('d',masses), array('d',values))
            getattr(ws, "import")(splineCentral, ROOT.RooFit.RecycleConflictNodes())
            shiftFormula = '{0}'.format(centralName)
            for shift in shifts:
                up = [u-c for u,c in zip(shifts[shift]['up'],values)]
                down = [d-c for d,c in zip(shifts[shift]['down'],values)]
                upName = '{0}_{1}Up'.format(splineName,shift)
                downName = '{0}_{1}Down'.format(splineName,shift)
                splineUp   = ROOT.RooSpline1D(upName,  upName,  ws.var(self.mh), len(masses), array('d',masses), array('d',up))
                splineDown = ROOT.RooSpline1D(downName,downName,ws.var(self.mh), len(masses), array('d',masses), array('d',down))
                getattr(ws, "import")(splineUp, ROOT.RooFit.RecycleConflictNodes())
                getattr(ws, "import")(splineDown, ROOT.RooFit.RecycleConflictNodes())
                shiftFormula += ' + TMath::Max(0,{shift})*{upName} + TMath::Min(0,{shift})*{downName}'.format(shift=shift,upName=upName,downName=downName)
            spline = ROOT.RooFormulaVar(splineName, splineName, shiftFormula, ROOT.RooArgList())
        else:
            spline = ROOT.RooSpline1D(splineName,  splineName,  ws.var(self.mh), len(masses), array('d',masses), array('d',values))
        getattr(ws, "import")(spline, ROOT.RooFit.RecycleConflictNodes())


class Gaussian(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        sigmaName = 'sigma_{0}'.format(label)
        mean = self.kwargs.get('mean',[1,0,1000])
        sigma = self.kwargs.get('sigma',[1,0,100])
        # variables
        ws.factory('{0}[{1}, {2}, {3}]'.format(meanName,*mean))
        ws.factory('{0}[{1}, {2}, {3}]'.format(sigmaName,*sigma))
        # build model
        ws.factory("Gaussian::{0}({1}, {2}, {3})".format(label,self.x,meanName,sigmaName))
        self.params = [meanName,sigmaName]

class GaussianSpline(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        sigmaName = 'sigma_{0}'.format(label)
        masses = self.kwargs.get('masses', [])
        means  = self.kwargs.get('means',  [])
        sigmas = self.kwargs.get('sigmas', [])
        # splines
        meanSpline  = ROOT.RooSpline1D(meanName,  meanName,  ws.var('MH'), len(masses), array('d',masses), array('d',means))
        sigmaSpline = ROOT.RooSpline1D(sigmaName, sigmaName, ws.var('MH'), len(masses), array('d',masses), array('d',sigmas))
        # import
        getattr(ws, "import")(meanSpline, ROOT.RooFit.RecycleConflictNodes())
        getattr(ws, "import")(sigmaSpline, ROOT.RooFit.RecycleConflictNodes())
        # build model
        ws.factory("Gaussian::{0}({1}, {2}, {3})".format(label,self.x,meanName,sigmaName))
        self.params = [meanName,sigmaName]

class BreitWigner(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        widthName = 'width_{0}'.format(label)
        mean  = self.kwargs.get('mean',  [1,0,1000])
        width = self.kwargs.get('width', [1,0,100])
        # variables
        ws.factory('{0}[{1}, {2}, {3}]'.format(meanName,*mean))
        ws.factory('{0}[{1}, {2}, {3}]'.format(widthName,*width))
        # build model
        ws.factory("BreitWigner::{0}({1}, {2}, {3})".format(label,self.x,meanName,widthName))
        self.params = [meanName,widthName]

class BreitWignerSpline(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        widthName = 'width_{0}'.format(label)
        masses = self.kwargs.get('masses', [])
        means  = self.kwargs.get('means',  [])
        widths = self.kwargs.get('widths', [])
        # splines
        meanSpline  = ROOT.RooSpline1D(meanName,  meanName,  ws.var('MH'), len(masses), array('d',masses), array('d',means))
        widthSpline = ROOT.RooSpline1D(widthName, widthName, ws.var('MH'), len(masses), array('d',masses), array('d',widths))
        # import
        getattr(ws, "import")(meanSpline, ROOT.RooFit.RecycleConflictNodes())
        getattr(ws, "import")(widthSpline, ROOT.RooFit.RecycleConflictNodes())
        # build model
        ws.factory("BreitWigner::{0}({1}, {2}, {3})".format(label,self.x,meanName,widthName))
        self.params = [meanName,widthName]

class Voigtian(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        widthName = 'width_{0}'.format(label)
        sigmaName = 'sigma_{0}'.format(label)
        mean  = self.kwargs.get('mean',  [1,0,1000])
        width = self.kwargs.get('width', [1,0,100])
        sigma = self.kwargs.get('sigma', [1,0,100])
        # variables
        ws.factory('{0}[{1}, {2}, {3}]'.format(meanName,*mean))
        ws.factory('{0}[{1}, {2}, {3}]'.format(widthName,*width))
        ws.factory('{0}[{1}, {2}, {3}]'.format(sigmaName,*sigma))
        # build model
        ws.factory("Voigtian::{0}({1}, {2}, {3}, {4})".format(label,self.x,meanName,widthName,sigmaName))
        self.params = [meanName,widthName,sigmaName]

class VoigtianSpline(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        widthName = 'width_{0}'.format(label)
        sigmaName = 'sigma_{0}'.format(label)
        masses = self.kwargs.get('masses', [])
        means  = self.kwargs.get('means',  [])
        widths = self.kwargs.get('widths', [])
        sigmas = self.kwargs.get('sigmas', [])
        # splines
        meanSpline  = ROOT.RooSpline1D(meanName,  meanName,  ws.var('MH'), len(masses), array('d',masses), array('d',means))
        widthSpline = ROOT.RooSpline1D(widthName, widthName, ws.var('MH'), len(masses), array('d',masses), array('d',widths))
        sigmaSpline = ROOT.RooSpline1D(sigmaName, sigmaName, ws.var('MH'), len(masses), array('d',masses), array('d',sigmas))
        # import
        getattr(ws, "import")(meanSpline, ROOT.RooFit.RecycleConflictNodes())
        getattr(ws, "import")(widthSpline, ROOT.RooFit.RecycleConflictNodes())
        getattr(ws, "import")(sigmaSpline, ROOT.RooFit.RecycleConflictNodes())
        # build model
        ws.factory("Voigtian::{0}({1}, {2}, {3}, {4})".format(label,self.x,meanName,widthName,sigmaName))
        self.params = [meanName,widthName,sigmaName]

class Exponential(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        lambdaName = 'lambda_{0}'.format(label)
        lamb = self.kwargs.get('lamb',  [-1,-5,0])
        # variables
        ws.factory('{0}[{1}, {2}, {3}]'.format(lambdaName,*lamb))
        # build model
        ws.factory("Exponential::{0}({1}, {2})".format(label,self.x,lambdaName))
        self.params = [lambdaName]

class Sum(Model):

    def __init__(self,name,**kwargs):
        self.doRecursive = kwargs.pop('recursive',False)
        self.doExtended = kwargs.pop('extended',False)
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        pdfs = []
        for n, (pdf, r) in enumerate(self.kwargs.iteritems()):
            ws.factory('{0}_norm[{1},{2}]'.format(pdf,*r))
            pdfs += [pdf]
        # build model
        if self.doRecursive:
            prev = ''
            sumargs = []
            for pdf in pdfs:
                curr = '{0}_norm*{1}*{0}'.format(pdf,prev) if prev else '{0}_norm*{0}'.format(pdf)
                sumargs += [curr]
                prev = '{0}*(1-{1}_norm)'.format(prev,pdf) if prev else '(1-{0}_norm)'.format(pdf)
        elif self.doExtended:
            sumargs = ['{0}_norm*{0}'.format(pdf) for pdf in pdfs]
        else: # Don't do this if you have more than 2 pdfs ...
            if len(pdfs)>2: logging.warning('This sum is not guaranteed to be positive because there are more than two arguments. Better to use the option recursive=True.')
            sumargs = ['{0}'.format(pdf) if len(pdfs)==n+1 else '{0}_norm*{0}'.format(pdf) for n,pdf in enumerate(pdfs)]
        ws.factory("SUM::{0}({1})".format(label, ', '.join(sumargs)))
        self.params = ['{}_norm'.format(pdf) for pdf in pdfs]
