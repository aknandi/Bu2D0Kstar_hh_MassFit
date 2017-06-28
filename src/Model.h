#ifndef Model_h
#define Model_h
#include "RooRealVar.h"
#include <string>

#include "RooSimultaneous.h"
#include "RooCategory.h"

#include "Base.h"
#include "Settings.h"
#include "Yields.h"
#include "Pdf_Fit.h"
#include "Pdf_Gen.h"
#include "RooDataSet.h"
#include "RooGaussian.h"

class Model : public Base
{
 public :
  Model(Settings*,RooRealVar*,RooCategory*,std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<std::string>);
  ~Model(){}
  Pdf_Fit fitPdf;
  Pdf_Gen genPdf;
  RooSimultaneous* getGenPdf();
  RooSimultaneous* getFitPdf();

  void printYieldsAndPurities(string, double, double, RooFitResult*);
  std::map< std::string,std::map< std::string, std::map< std::string,std::map< std::string, std::map<std::string,double > > > > > plotNums;

  std::vector <RooRealVar*> *GetFixedParameters() {return fitPdf.GetFixedParameters();}

  double totalSignificance = 0;
  double errSignificance = 0;

 private:
  Yields* yields;

  Settings* _genConfs;
  RooSimultaneous* sim;
	
  RooRealVar* mB;
  RooCategory* _cat;
	
  RooGaussian* gaus_f010;
  RooGaussian* gaus_r_dstkst;

  std::vector<std::string> _modeList;
  std::vector<std::string> _chargeList;
  std::vector<std::string> _trackList;
  std::vector<std::string> _runList;
  
};


#endif
