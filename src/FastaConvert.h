#ifndef FASTA_CONVERT
#define FASTA_CONVERT

#include <map>
#include <string>

class FastaConvert {
 public :
  FastaConvert();
  
  void execute(std::map<std::string, std::string> & parameters);
};

#endif //FASTA_CONVERT
