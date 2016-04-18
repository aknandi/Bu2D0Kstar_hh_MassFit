#include "Settings.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "boost/algorithm/string.hpp"

std::vector<std::pair<std::string,std::string> > Settings::readPairStringsToVector(std::string filename)
{
  bool db(false);
  std::string line;	// string to read into from file
  std::ifstream myfile;
 
  std::vector<std::pair<std::string,std::string> > return_vec; // Vector to return
  // Now try and open the file
  try {
    myfile.open(filename.c_str());
    while(myfile.good() && !myfile.eof()) {
      getline(myfile,line);
      boost::trim(line);//Trim away any whitespace (so that it can't be interpreted as 'zero')
      if(myfile.eof())	break; // If at the end of file, break out
      if(line[0]=='*') continue; // If line is commented, skip it
      if(line=="") continue; // If line is empty, skip it
      // Now got the line, need to split it (whether it's a tab or a space)
      int pos_sp = line.find(' ');
      int pos_tab = line.find('\t');
      if(pos_sp==(int)std::string::npos) pos_sp=pos_tab;
      std::pair<std::string,std::string> tempPair;
      // Should have at least two entries in the string
      if(pos_sp==(int)std::string::npos) {
        std::cout << "Parameter file " << filename << " needs at least two items per line for readPairStringsToVector() to work" << std::endl;
        std::exit(2);
      }
      else {
        tempPair.first=line.substr(0,pos_sp);
        tempPair.second=line.substr(pos_sp+1,line.length()-pos_sp-1);
        // Now trim away any whitespace
        boost::trim(tempPair.first);
        boost::trim(tempPair.second);
      }
      return_vec.push_back(tempPair);
    }
  }
  catch (std::ifstream::failure e) {
    std::cout << "Exception opening/reading file" << std::endl;
  }

  if(db) {
    for (std::vector<std::pair<std::string,std::string> >::iterator it=return_vec.begin(); it!=return_vec.end(); it++) {
      //			std::cout << "Element 1: " << it->first << " and 2: " << it->second << std::endl;
    }
  }
  return return_vec;
}

std::map<std::string,std::string> Settings::readPairStringsToMap(std::string filename)
{
  // Use the function above and return a map
  std::vector<std::pair<std::string,std::string> > settings_vec = readPairStringsToVector(filename); 
  std::map<std::string,std::string> settings_map;
  for (std::vector<std::pair<std::string,std::string> >::iterator it=settings_vec.begin(); it!=settings_vec.end(); it++) {
    if(settings_map.find(it->first)!=settings_map.end()) {
      std::cout << "Entry (" << it->first << ","<<it->second << ") in file " << filename << " already in map. Exiting - go figure" << std::endl;
      std::exit(2);
    }
    else {
      settings_map[it->first]=it->second;
      //debug line std::cout<<it->first<< " "<<it->second<<std::endl;		

    }
  }
  _pairStringsMap.insert(settings_map.begin(),settings_map.end());

  _haveReadInMap=true;
  return settings_map;
}

std::string Settings::get(std::string requestedKey)
{
  if(_haveReadInMap==false) {
    std::cout << "Need to read into a map before calling the special map get function, so exiting. " << std::endl;
    exit(2);
  }
  if(this->contains(requestedKey)==false) {
    std::cout << "Settings object with name \"" << _myName << "\" does not contain key: " << requestedKey << ". Exiting!" <<std::endl;
    exit(2);
  }
  return _pairStringsMap[requestedKey];
}

double Settings::getD(std::string requestedKey)
{
  //std::cout << "Returning from getD a value for "<<requestedKey<< " of "<<atof(get(requestedKey).c_str()) << std::endl;
  return atof(get(requestedKey).c_str());
}

int Settings::getI(std::string requestedKey)
{
  return atoi(get(requestedKey).c_str());
}

bool Settings::isChargeSeparated() {
	return get("chargeSeparated") == "true";
}
 
bool Settings::contains(std::string requestedKey)
{
  if(_haveReadInMap==false) {
    std::cout << "Need to read into a map before calling the special map contains function, so exiting. " << std::endl;
    exit(2);
  }
  if (_pairStringsMap.find(requestedKey)==_pairStringsMap.end()) return false;
  return true;
}

