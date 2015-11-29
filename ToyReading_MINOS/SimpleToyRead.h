#ifndef SimpleToyRead_h
#define SimpleToyRead_h

#include <string>
#include "TRandom3.h"
#include "TFile.h"
#include <iostream>

class SimpleToyRead {

 public :
  SimpleToyRead();
 

  void MakeNtupleFromTextFile(std::string name_var, float trueval);
  void MakeSomePlotsFromRootFile(std::string name_var, float trueval);

  std::ofstream outputfile;
};

#endif

