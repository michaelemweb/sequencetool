#include "Utils.h"

#include <iostream>
#include <stdlib.h>

bool require(const std::vector<Parameter> & keys, 
	     const std::map<std::string, std::string> & parameters)
{
  bool returnVal = true;

  unsigned i;
  for (i = 0; i < keys.size(); i++) {
    std::string key = keys[i].name_;
    std::string description = keys[i].description_;
    
    std::map<std::string, std::string>::const_iterator i = parameters.find(key);

    if (i == parameters.end()) {
      std::cerr << "Missing module argument: " << key << std::endl;
      if (description != "")
	std::cerr << description << std::endl;
      returnVal = false;
    }
  }

  return returnVal;
}


void handleParseException(seq::ParseException& e)
{
  if (e.recovered())
    std::cout << e.name() << ",\"" << e.message() << "\",,,,,"
	      << std::endl;
  else
    std::cerr << "Fatal error: " << e.message() << std::endl;
  
  if (!e.recovered())
    exit(1);
}
