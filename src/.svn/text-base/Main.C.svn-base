#include <iostream>
#include <string>

#include "TApplication.h"
#include "TROOT.h"
#include "RooMsgService.h"

#include "Settings.h"
#include "Fitting.h"

#include "boost/algorithm/string.hpp"

int main(int argc, char *argv[])
{
  if(argc!=2)
    {
      cout << "Please provide a filename like Settings/GeneralSettings.txt" << endl;
      exit(2);
    }
  gROOT->SetStyle("Plain");
	
  TApplication app("app",0,0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  Settings* generalSettings = new Settings();
  generalSettings->readPairStringsToMap(argv[1]);

  // Check sanity of settings:
  std::string batchMode = generalSettings->get("batchMode");boost::to_lower(batchMode);
  if(batchMode=="true"){cout << "OPERATING IN BATCH MODE" << endl;gROOT->SetBatch(kTRUE);}
  Fitting fit(&app,generalSettings);
}
