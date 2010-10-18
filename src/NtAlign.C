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

void writeSequenceWithoutInsertions(std::ostream& o,
				    const seq::NTSequence& reference,
				    const seq::NTSequence& target)
{
  o << ">" << target.name() << " " << target.description() << std::endl;

  for (unsigned i = 0; i < reference.size(); ++i)
    if (reference[i] != seq::Nucleotide::GAP)
      o << target[i];

  o << std::endl;
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

  const int cutoff = atoi(parameters["cutoff"].c_str());

  std::ifstream reference_f(parameters["reference"].c_str());
  std::ifstream target_f(parameters["target"].c_str());

  std::ofstream output_f(parameters["output"].c_str());

  try {
    seq::NTSequence reference;
    reference_f >> reference;
    
    while (target_f) {
      seq::NTSequence s;
      target_f >> s;
      if (target_f) {
	seq::NTSequence r = reference;
	double score = seq::Align(r, s);
	if (score >= cutoff)
	  writeSequenceWithoutInsertions(output_f, r, s);
      }
    }
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

