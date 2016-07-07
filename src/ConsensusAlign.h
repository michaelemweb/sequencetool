#ifndef CONSENSUS_ALIGN
#define CONSENSUS_ALIGN

#include <map>
#include <string>

class ConsensusAlign {
public :
  ConsensusAlign();
  
  void execute(std::map<std::string, std::string> & parameters);
};

#endif //CONSENSUS_ALIGN
