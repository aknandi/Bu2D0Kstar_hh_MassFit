#include "Settings.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "InternalStorage.h"
//#include "CommonTools.h"
using namespace std;

void InternalStorage::Initialize(Settings* Input , std::vector<std::string> allmodeList, std::vector<std::string> allchargeList, std::vector<std::string> alltrackList, std::vector<std::string> allbinList)
{

  _Control=Input;

  modeList=allmodeList;
  chargeList = allchargeList;
  trackList=alltrackList;
  binList=allbinList;

  
  SetTotals();
  std::cout<<" I've set the Totals "<<std::endl;

}

InternalStorage::InternalStorage()
{

}


void InternalStorage::SetTotals()
{
  std::string filename = _Control->get("PathnameToTotals");

  Settings totals;
  totals.readPairStringsToMap(filename);
 
  for(std::vector<std::string>::iterator m=modeList.begin(); m!=modeList.end(); m++){
    for(std::vector<std::string>::iterator c=chargeList.begin(); c!=chargeList.end(); c++){
      for(std::vector<std::string>::iterator t=trackList.begin(); t!=trackList.end(); t++){
        /*
        FracBsInLowMass[*m]["pass"][*c][*t]=totals.getD("FracBsBkgInLowMass_"+(*m)+"_pass_"+(*c)+"_"+(*t));
          std::string SigYield = "Nsignal_"+(*m)+"_"+(*p)+"_"+(*c)+"_"+(*t);
          std::string ComYield = "NComb_"+(*m)+"_"+(*p)+"_"+(*c)+"_"+(*t);
          std::string LowYield = "NLow_"+(*m)+"_"+(*p)+"_"+(*c)+"_"+(*t);
          if(!totals.contains(SigYield)){std::cout<<"Missing key: "<<SigYield<<std::endl;exit(2);}
          if(!totals.contains(ComYield)){std::cout<<"Missing key: "<<ComYield<<std::endl;exit(2);}
          if(!totals.contains(LowYield)){std::cout<<"Missing key: "<<LowYield<<std::endl;exit(2);}

          TotSig[*m][*p][*c][*t]=totals.getD(SigYield);
          TotLowmass[*m][*p][*c][*t]=totals.getD(LowYield);
          TotComb[*m][*p][*c][*t]=totals.getD(ComYield);
          if(totals.getI("DEBUG_runInHighStatsMode")) {
            TotSig[*m][*p][*c][*t]*=100.;
            TotLowmass[*m][*p][*c][*t]*=100.;
            TotComb[*m][*p][*c][*t]*=100.;
          }
          */
      }
    }
  }
}

   
