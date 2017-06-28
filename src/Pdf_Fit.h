#ifndef Pdf_Fit_h
#define Pdf_Fit_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "Settings.h"
#include "Pdf_Base.h"
#include "myGaussian.h"
#include "DoubleCrystalBall.h"
#include "Exponential.h"
#include "PartRecoDstKst.h"
#include "RooKeysPdf.h"
#include "myCruijff.h"

class Pdf_Fit
{
 public :
  Pdf_Fit(){};
  Pdf_Fit(Settings*,Settings*,RooRealVar*, std::vector<std::string>, std::vector<std::string>, std::vector<std::string> , std::vector<std::string>, int systematicFactor=0, std::string MCsimfit="");
  void setRelations();
  // PDFs to return
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bu;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_comb;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_lckst;

  std::vector <RooRealVar*> *GetFixedParameters() {return fixedParams;}

 private:
  Settings* _fileList;
  Settings* _genConfs;
  std::string _MCsimfit;
  std::string _toys;
  // PDFs
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,DoubleCrystalBall*> > > > bu;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Exponential*> > > > comb;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,PartRecoDstKst*> > > > dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,myCruijff*> > > > lckst;

  std::vector<std::string> _modeList;
  std::vector<std::string> _chargeList;
  std::vector<std::string> _trackTypeList;
  std::vector<std::string> _runList;

  //fixed parameters
  std::vector <RooRealVar*> *fixedParams;
};

#endif

