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

class BreitWignerSpline(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        meanName = 'mean_{0}'.format(label)
        sigmaName = 'sigma_{0}'.format(label)
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

class Sum(Model):

    def __init__(self,name,**kwargs):
        super(self.__class__,self).__init__(name,**kwargs)

    def build(self,ws,label):
        pdfs = []
        for n, (pdf, r) in enumerate(self.kwargs.iteritems()):
            ws.factory('{0}n{1}[{2},{3}]'.format(label,n,*r))
            pdfs += [pdf]
        # build model
        ws.factory("SUM::{0}({1})".format(label, ', '.join([pdf if len(pdf)==n else '{0}n{1}*{2}'.format(label,n,pdf) for n,pdf in enumerate(pdfs)])))
