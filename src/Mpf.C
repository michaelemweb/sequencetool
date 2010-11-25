#include <map>
#include <string>

#include <iostream>
#include <stdlib.h>

#include "NtAlign.h"
#include "FastaConvert.h"

int main(int argc, char **argv) {
  unsigned int i;
  
  int obligatoryParams = 1;
  if(argc < obligatoryParams + 1) {
    std::cerr << "Please provide a module." << std::endl;
    std::cerr << "The following modules are available: "
	      << "nt_align "
	      << "fasta-convert " << std::endl;
    exit(0);
  }
  std::string module = argv[1];
  
  int amountOfParameters = argc - obligatoryParams - 1;
  if (amountOfParameters%2 == 1) {
    std::cerr << "Please provide module parameters as:"
	      << " --parameterName parameterValue" 
	      << std::endl;	
    exit(0);
  }
  
  std::map<std::string, std::string> parameters;
  for(i = obligatoryParams+1; 
      i < amountOfParameters + obligatoryParams; 
      i=i+2) {
    std::string key = argv[i];
    std::string value = argv[i + 1];
    if (key.substr(0, 2) == "--")
      parameters[key.substr(2)] = value;
    else {
      std::cerr << "A module parameters should start with --" << std::endl;
      exit(0);
    }
  }

  if (module == "nt_align") {
    NtAlign align;
    align.execute(parameters);
  } else if (module == "fasta-convert") {
    FastaConvert convert;
    convert.execute(parameters);
  }
}
