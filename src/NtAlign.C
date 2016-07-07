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

#include "ssw_cpp.h"

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

void getQueryWithoutInsertions(seq::NTSequence& ref, seq::NTSequence& query,
			       int ref_begin,
			       int ref_end,
			       int query_begin,
			       int query_end,
			       const std::vector<uint32_t>& cigar)
{
  /* truncate query to part that was aligned */
  if (query_end < query.size())
    query.erase(query.begin() + query_end, query.end());
    
  if (query_begin > 0)
    query.erase(query.begin(), query.begin() + query_begin);

  unsigned pos = 0;
  for (unsigned i = 0; i < cigar.size(); ++i) {
    int length = cigar[i] >> 4;
    int op = cigar[i] & 0xF;

    switch (op) {
    case 0: // M
    case 7: // =
    case 8: // X
      pos += length;
      break;
    case 1: // I
      query.erase(query.begin() + pos, query.begin() + pos + length);
      break;
    case 2: // D
    case 3: // N
      for (unsigned j = 0; j < length; ++j)
	query.insert(query.begin() + pos, seq::Nucleotide::GAP);
      pos += length;
      break;
    default:
    case 4: // S
      // ignoring since that seems to be quite inconsistent?
      break;
    }
  }

  /* insert leading and ending '-' */
  for (unsigned i = 0; i < ref_begin; ++i)
    query.insert(query.begin(), seq::Nucleotide::GAP);

  for (unsigned i = ref_end; i < ref.size(); ++i)
    query.push_back(seq::Nucleotide::GAP);
}

double sswAlign(seq::NTSequence& ref, seq::NTSequence& target)
{
  std::string r = ref.asString();
  std::string t = target.asString();

  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  aligner.Align(t.c_str(), r.c_str(), r.size(), filter, &alignment);

  getQueryWithoutInsertions(ref, target, alignment.ref_begin,
			    alignment.ref_end,
			    alignment.query_begin,
			    alignment.query_end,
			    alignment.cigar);

  return alignment.sw_score;
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
	int end = reference.size();

	/*
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
	*/

	seq::NTSequence r(reference.begin() + begin, reference.begin() + end);

	if (r.size() * s.size() < 5000 * 5000) {
	  seq::NTSequence r1 = r;
	  seq::NTSequence s1 = s;
	  double score1 = seq::Align(r1, s1);
	  seq::NTSequence r2 = r;
	  seq::NTSequence s2 = s.reverseComplement();
	  double score2 = seq::Align(r2, s2);
	  if (std::max(score1, score2) >= cutoff) {
	    bool use1 = score1 > score2;
	    seq::NTSequence& rr = use1 ? r1 : r2;
	    seq::NTSequence& ss = use1 ? s1 : s2;
	    writeSequenceWithoutInsertions(output_f, rr, ss, 
					   begin, reference.size());
	    if (use1)
	      std::cerr << s;
	    else
	      std::cerr << s.reverseComplement();
	  }
	} else {
	  seq::NTSequence r1 = r;
	  seq::NTSequence s1 = s;
	  double score1 = sswAlign(r1, s1);
	  seq::NTSequence r2 = r;
	  seq::NTSequence s2 = s.reverseComplement();
	  double score2 = sswAlign(r2, s2);
	  if (std::max(score1, score2) >= cutoff) {
	    bool use1 = score1 > score2;
	    seq::NTSequence& ss = use1 ? s1 : s2;
	    output_f << ss;

	    if (use1)
	      std::cerr << s;
	    else
	      std::cerr << s.reverseComplement();
	  }
	}
      }
    }
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

