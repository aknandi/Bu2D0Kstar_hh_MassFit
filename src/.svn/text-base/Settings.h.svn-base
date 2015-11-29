#ifndef Settings_h
#define Settings_h

#include <string>
#include <vector>
#include <map>
#include <utility>

class Settings
{
 public:
  Settings(){_haveReadInMap=false;_myName="Settings::NameNotSet";}
  Settings(std::string myName){_haveReadInMap=false;_myName=myName;}

  std::vector< std::pair<std::string,std::string> > readPairStringsToVector(std::string);
  std::map< std::string, std::string> readPairStringsToMap(std::string);
		
  // And a special 'get' function to save a lot of typing with accessing const maps ("myMap.find("key")->second" => "myMap.get("key"))
  std::string get(std::string);
  double getD(std::string); // uses 'get' and returns a double
  int getI(std::string); // uses 'get' and returns an integer

  // check if key contained before query map
  bool contains(std::string);
 private:
  std::map< std::string, std::string> _pairStringsMap;
  bool _haveReadInMap;
  std::string _myName; // useful for error reporting ('Settings object with name xxx didn't contain requested key')
};

#endif
