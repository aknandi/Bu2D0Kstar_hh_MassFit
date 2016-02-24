#ifndef Pdf_Fit_h
#define Pdf_Fit_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "Settings.h"
#include "Pdf_Base.h"
#include "myGaussian.h"
#include "myCrystalBall.h"
#include "Exponential.h"
#include "PartRecoDstKst.h"

class Pdf_Fit
{
 public :
  Pdf_Fit();
  Pdf_Fit(Settings*,Settings*,RooRealVar*, std::vector<std::string>, std::vector<std::string>, std::vector<std::string> , std::vector<std::string>, int systematicFactor=0, std::string MCsimfit="");
  void setRelations();
  // PDFs to return
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bu;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bs;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_comb;
/*  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_drho;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bd_dstkst;*/
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bu_dstkst;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bs_dstkst_001;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bs_dstkst_010;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_dkpipi;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_dpipipi;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_lambda;
//  std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > gaus_frac010_bs;
//  std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > gaus_frac010_bd;
  //RooAbsPdf* gaus_frac010_bs;
  //RooAbsPdf* gaus_frac010_bd;
  //RooRealVar* bs_frac010;

  std::vector <RooRealVar*> *GetFixedParameters() {return fixedParams;}

 private:
  Settings* _fileList;
  Settings* _genConfs;
  std::string _MCsimfit;
  std::string _toys;
  // PDFs
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,myCrystalBall*> > > > bu;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bs;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Exponential*> > > > comb;
//  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > drho;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,PartRecoDstKst*> > > > bu_dstkst;
/*  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bs_dstkst_001;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bs_dstkst_010;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bd_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > dkpipi;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > dpipipi;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > lambda;*/

  std::vector<std::string> _modeList;
  std::vector<std::string> _chargeList;
  std::vector<std::string> _trackTypeList;
  std::vector<std::string> _runList;

  //fixed parameters
  std::vector <RooRealVar*> *fixedParams;
};

#endif

