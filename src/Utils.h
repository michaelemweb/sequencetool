#ifndef UTILS
#define UTILS

#include <vector>
#include <map>
#include <string>

#include "ParseException.h"

class Parameter {
public:
 Parameter(std::string name) : 
  name_(name)
  {}

public:
  std::string name_;
  std::string description_;
};

bool require(const std::vector<Parameter> & keys, 
	     const std::map<std::string, std::string> & parameters);

void handleParseException(seq::ParseException& e);

#endif //UTILS
