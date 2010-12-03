#include "NtAlign.h"

#include "NTSequence.h"
#include "Align.h"

#include "Utils.h"

#include <boost/algorithm/string.hpp>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdlib.h>

NtAlign::NtAlign()
{
}

void writeSequenceWithoutInsertions(std::ostream& o,
				    const seq::NTSequence& reference,
				    const seq::NTSequence& target,
				    int begin,
				    int size)
{
  o << ">" << target.name() << " " << target.description() << std::endl;

  int j = 0;
  for (unsigned i = 0; i < begin; ++i) {
    o << "-";
    ++j;
  }

  for (unsigned i = 0; i < reference.size(); ++i)
    if (reference[i] != seq::Nucleotide::GAP) {
      o << target[i];
      ++j;
    }

  for (unsigned i = 0; i < (size - j); ++i)
    o << "-";

  o << std::endl;
}

void NtAlign::execute(std::map<std::string, std::string> & parameters)
{
  std::vector<Parameter> requiredKeys;
  requiredKeys.push_back(Parameter("reference"));
  requiredKeys.push_back(Parameter("target"));
  requiredKeys.push_back(Parameter("output"));
  requiredKeys.push_back(Parameter("cutoff"));
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
	int begin = 0;
	int end = reference.size() - 1;

	std::string d = s.description();
	boost::trim(d);
	if (!d.empty()) {
	  std::vector<std::string> values;
	  boost::split(values, d, boost::is_any_of("-"));

	  if (values.size() == 2) {
	    const int MARGIN = 50;
	    begin = std::max(begin,
			     atoi(values[0].substr(1).c_str()) - MARGIN);
	    end = std::min(end,
			   atoi(values[1].substr(0, values[1].length() - 1)
				.c_str()) + MARGIN);
	  }
	}

	seq::NTSequence r(reference.begin() + begin, reference.begin() + end);

	double score = seq::Align(r, s);
	if (score >= cutoff)
	  writeSequenceWithoutInsertions(output_f, r, s, begin, reference.size());
      }
    }
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

