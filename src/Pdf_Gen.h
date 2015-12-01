#ifndef Pdf_Gen_h
#define Pdf_Gen_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "Settings.h"
#include "Pdf_Base.h"

class Pdf_Gen
{
 public :
  Pdf_Gen(){};
  Pdf_Gen(Settings*,RooRealVar*, std::vector<std::string>, std::vector<std::string>, std::vector<std::string> , std::vector<std::string>);
  void setRelations();
  // PDFs to return
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bu;
/*  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bs;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_comb;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_drho;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bd_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,RooAbsPdf*> > > > roopdf_bs_dstkst;*/

  std::vector <RooRealVar*> *GetFixedParameters() {return fixedParams;}

 private:
  Settings* _fileList;
  std::string _MCsimfit;
  // PDFs
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bu;
/*  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bs;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > comb;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > drho;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bs_dstkst;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string,Pdf_Base*> > > > bd_dstkst;*/

  std::vector<std::string> _modeList;
  std::vector<std::string> _chargeList;
  std::vector<std::string> _trackTypeList;
  std::vector<std::string> _runList;

  //fixed parameters
  std::vector <RooRealVar*> *fixedParams;
};

#endif

