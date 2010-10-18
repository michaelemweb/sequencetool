#ifndef UTILS
#define UTILS

#include <set>
#include <map>
#include <string>
#include <iostream>

bool require(const std::set<std::string> & keys, 
	     const std::map<std::string, std::string> & parameters)
{
  std::set<std::string>::iterator key_it;

  bool returnVal = true;

  for (key_it = keys.begin(); key_it != keys.end(); key_it++) {
    std::string key = *key_it;
    std::map<std::string, std::string>::const_iterator i = parameters.find(key);

    if (i == parameters.end()) {
      std::cerr << "Missing module argument: " << key << std::endl;
      returnVal = false;
    }
  }

  return returnVal;
}

#endif //UTILS
