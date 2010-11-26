#include "FastaConvert.h"

#include "NTSequence.h"

#include "Utils.h"

#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdlib.h>

FastaConvert::FastaConvert()
{
}

void FastaConvert::execute(std::map<std::string, std::string> & parameters)
{
  std::vector<Parameter> requiredKeys;
  requiredKeys.push_back(Parameter("input-fasta"));
  requiredKeys.push_back(Parameter("output-file"));
  Parameter outputFormat("output-format");
  outputFormat.description_ = "Possible formats are: stockholm.";
  requiredKeys.push_back(outputFormat);
  if (!require(requiredKeys, parameters))
    exit(1);

  std::ifstream input_f(parameters["input-fasta"].c_str());

  std::ofstream output_f(parameters["output-file"].c_str());
  std::string format = parameters["output-format"];

  std::set<seq::NTSequence> sequences;
  try {
    while (input_f) {
      seq::NTSequence s;
      input_f >> s;
      if (input_f)
	sequences.insert(s);
    }
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }

  if (format == "stockholm") {
    const int lineLenght = 10000;
    seq::writeStockholm(output_f, sequences, lineLenght);
  } else {
    std::cerr << format << " is not a supported output format" << std::endl;
    exit(0);
  }
}

