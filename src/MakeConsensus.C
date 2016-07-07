#include "MakeConsensus.h"

#include "NTSequence.h"
#include "Utils.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdlib.h>

MakeConsensus::MakeConsensus()
{
}

struct Contig {
  seq::NTSequence seq;
  double cov;
  int start, end;
};

int firstNonGap(const seq::NTSequence& seq)
{
  for (unsigned i = 0; i < seq.size(); ++i)
    if (seq[i] != seq::Nucleotide::GAP)
      return i;
  return -1;
}

int lastNonGap(const seq::NTSequence& seq)
{
  for (int i = seq.size() - 1; i >= 0; --i)
    if (seq[i] != seq::Nucleotide::GAP)
      return i;
  return -1;
}

void writeContig(std::ostream& o, seq::NTSequence& contig, int index)
{
  contig.setName("contig_" + boost::lexical_cast<std::string>(index)
		 + "_len_" + boost::lexical_cast<std::string>(contig.size()));
  o << contig;
}

void MakeConsensus::execute(std::map<std::string, std::string> & parameters)
{
  std::vector<Parameter> requiredKeys;
  requiredKeys.push_back(Parameter("input"));
  requiredKeys.push_back(Parameter("output"));
  requiredKeys.push_back(Parameter("max-gap"));
  requiredKeys.push_back(Parameter("max-missing"));
  requiredKeys.push_back(Parameter("min-count"));
  if (!require(requiredKeys, parameters))
    exit(1);

  const int maxGap = atoi(parameters["max-gap"].c_str());
  const int maxMissing = atoi(parameters["max-missing"].c_str());
  const double minCount = atof(parameters["min-count"].c_str());

  std::ifstream input_f(parameters["input"].c_str());
  std::ofstream output_f(parameters["output"].c_str());

  std::vector<Contig> alignment;

  try {
    while (input_f) {
      seq::NTSequence s;
      input_f >> s;
      if (input_f) {
	Contig c;
	c.seq = s;
	c.cov = 1;
	std::size_t covPos = s.name().find("_cov_");
	if (covPos != std::string::npos)
	  c.cov = atof(s.name().substr(covPos + 5).c_str());
	c.start = firstNonGap(s);
	c.end = lastNonGap(s);
	alignment.push_back(c);
      }
    }

    int contigs = 0;
    seq::NTSequence contig;
    const seq::NTSequence& ref = alignment[0].seq;

    int gapSize = 0;
    for (unsigned i = 0; i < ref.size(); ++i) {
      double base_counts[4]; // acgt
      double cov = 0;

      for (int l = 0; l < 4; ++l)
	base_counts[l] = 0.0;

      for (unsigned j = 1; j < alignment.size(); ++j) {
	const Contig& c = alignment[j];
	if (i >= c.start && i <= c.end) {
	  std::vector<seq::Nucleotide> nucs;
	  c.seq[i].nonAmbiguousNucleotides(nucs);
	  for (unsigned k = 0; k < nucs.size(); ++k) {
	    if (nucs[k].intRep() < 4)
	      base_counts[nucs[k].intRep()] += c.cov;
	  }

	  cov += c.cov;
	}
      }

      if (base_counts[0] + base_counts[1] + base_counts[2] + base_counts[3]
	  == 0) {
	if (!contig.empty()) {
	  ++gapSize;
	  if (gapSize > maxMissing && cov < minCount) {
	    writeContig(output_f, contig, ++contigs);
	    contig = seq::NTSequence();
	  }
	}
      } else {
	if (gapSize <= maxGap)
	  gapSize = 0;
	else {
	  for (int i = 0; i < gapSize; ++i)
	    contig.push_back(seq::Nucleotide::N);
	  gapSize = 0;
	}

	std::set<seq::Nucleotide> have, have_0;

	for (int i = 0; i < 4; ++i) {
	  if (base_counts[i] >= minCount)
	    have.insert(seq::Nucleotide::fromRep(i));
	  if (base_counts[i] > 0)
	    have_0.insert(seq::Nucleotide::fromRep(i));
	}

	if (have.empty())
	  have = have_0;

	seq::Nucleotide result = seq::Nucleotide::singleNucleotide(have);
	contig.push_back(result);
      }
    }

    if (!contig.empty())
      writeContig(output_f, contig, ++contigs);
  } catch(seq::ParseException& e) {
    handleParseException(e);
  }
}

