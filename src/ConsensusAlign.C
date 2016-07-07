#include "ConsensusAlign.h"

#include "NTSequence.h"
#include "NeedlemanWunsh.h"

#include "Utils.h"

#include <boost/algorithm/string.hpp>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdlib.h>

#include "ssw_cpp.h"

ConsensusAlign::ConsensusAlign()
{
}

void cigarToAlignment(seq::NTSequence& ref, seq::NTSequence& query,
		      int ref_begin,
		      int ref_end,
		      int query_begin,
		      int query_end,
		      const std::vector<uint32_t>& cigar)
{
  std::cerr << "ref_begin: " << ref_begin << std::endl;
  std::cerr << "ref_end: " << ref_end << std::endl;
  std::cerr << "query_begin: " << query_begin << std::endl;
  std::cerr << "query_end: " << query_end << std::endl;

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
      for (unsigned j = 0; j < length; ++j)
	ref.insert(ref.begin() + ref_begin + pos, seq::Nucleotide::GAP);
      pos += length;
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

  while (query.size() < ref.size())
    query.push_back(seq::Nucleotide::GAP);
}

double sswAlignFull(seq::NTSequence& ref, seq::NTSequence& target)
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

  cigarToAlignment(ref, target, alignment.ref_begin,
		   alignment.ref_end,
		   alignment.query_begin,
		   alignment.query_end,
		   alignment.cigar);

  return alignment.sw_score;
}

void mergeAlign(seq::NTSequence& r, seq::NTSequence& s,
		std::vector<seq::NTSequence>& alignment)
{
  if (alignment.empty()) {
    alignment.push_back(r);
    alignment.push_back(s);
  } else {
    seq::NTSequence& r1 = alignment[0];
    seq::NTSequence& r2 = r;
    unsigned i2 = 0;
    for (unsigned i1 = 0; i1 < r1.size();) {
      /* We need to do something special if gaps are different in r vs rr */
      if ((r1[i1] == seq::Nucleotide::GAP) !=
	  (r2[i2] == seq::Nucleotide::GAP)) {
	if (r1[i1] == seq::Nucleotide::GAP) {
	  // '-' vs 'n'
	  s.insert(s.begin() + i1, seq::Nucleotide::GAP);
	  ++i1;
	} else {
	  // 'n' vs '-'
	  for (unsigned j = 0; j < alignment.size(); ++j) {
	    seq::NTSequence& sa = alignment[j];
	    sa.insert(sa.begin() + i1, seq::Nucleotide::GAP);
	  }
	  ++i1;
	  ++i2;
	}
      } else {
	++i1;
	++i2;
      }
    }

    alignment.push_back(s);
  }
}

void ConsensusAlign::execute(std::map<std::string, std::string> & parameters)
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

  std::vector<seq::NTSequence> alignment;
  seq::NeedlemanWunsh nmw;

  try {
    seq::NTSequence reference;
    reference_f >> reference;

    while (target_f) {
      seq::NTSequence s;
      target_f >> s;
      if (target_f) {
	int begin = 0;
	int end = reference.size();

	seq::NTSequence r(reference.begin() + begin, reference.begin() + end);
	seq::NTSequence r1 = r;
	seq::NTSequence s1 = s;
	seq::NTSequence r2 = r;
	seq::NTSequence s2 = s.reverseComplement();

	double score1, score2;
	if (false && (r.size() * s.size() < 5000 * 5000)) {
	  score1 = nmw.align(r1, s1);
	  score2 = nmw.align(r2, s2);
	} else {
	  score1 = sswAlignFull(r1, s1);
	  score2 = sswAlignFull(r2, s2);
	}

	double score = std::max(score1, score2);
	if (score >= cutoff) {
	  std::cerr << s2.name() << ": " << score << std::endl;
	  bool use1 = score1 > score2;
	  seq::NTSequence& rr = use1 ? r1 : r2;
	  seq::NTSequence& ss = use1 ? s1 : s2;

	  // std::cout << rr << ss;

	  mergeAlign(rr, ss, alignment);
	}
      }
    }
    for (unsigned i = 0; i < alignment.size(); ++i)
      output_f << alignment[i];
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

