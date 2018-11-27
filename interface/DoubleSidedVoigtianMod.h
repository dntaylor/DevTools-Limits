#ifndef MY_DOUBLE_Voigtian
#define MY_DOUBLE_Voigtian

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooGlobalFunc.h"

#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"

#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooPlot.h"
 
class DoubleSidedVoigtianMod : public RooAbsPdf {
public:
  DoubleSidedVoigtianMod() {} ; 
  DoubleSidedVoigtianMod(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mean,
	      RooAbsReal& _sig1,
	      RooAbsReal& _sig2,
              RooAbsReal& _wid1,
              RooAbsReal& _wid2,
              Double_t _yMax);
  DoubleSidedVoigtianMod(const DoubleSidedVoigtianMod& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new DoubleSidedVoigtianMod(*this,newname); }
  inline virtual ~DoubleSidedVoigtianMod() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy sig1 ;
  RooRealProxy sig2 ;
  RooRealProxy wid1 ;
  RooRealProxy wid2 ;
  Double_t evaluate() const ;

private:

  Double_t yMax;
  ClassDef(DoubleSidedVoigtianMod,1) // Your description goes here...
};
 
#endif

