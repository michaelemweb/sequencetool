#include "NtAlign.h"

#include "NTSequence.h"
#include "Align.h"

#include "Utils.h"

#include <set>
#include <iostream>
#include <fstream>
#include <stdlib.h>

NtAlign::NtAlign()
{
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

void NtAlign::execute(std::map<std::string, std::string> & parameters)
{
  std::set<std::string> requiredKeys;
  requiredKeys.insert("reference");
  requiredKeys.insert("target");
  requiredKeys.insert("output");
  requiredKeys.insert("cutoff");
  if (!require(requiredKeys, parameters))
    exit(1);

  std::cerr << "reference" << parameters["reference"].c_str() << std::endl;
  std::cerr << "target" << parameters["target"].c_str() << std::endl;

  std::ifstream reference_f(parameters["reference"].c_str());
  std::ifstream target_f(parameters["target"].c_str());

  try {
    seq::NTSequence reference;
    reference_f >> reference;
    std::cout << reference;
    
    while (target_f) {
      seq::NTSequence s;
      target_f >> s;
      if (target_f) {
	//TODO
	//need to copy ref seq, otherwise after aligning with a target with insertions, the seq will contain gaps 
	double d = seq::Align(reference, s);
	
      }
    }
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

