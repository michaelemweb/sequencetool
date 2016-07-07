#ifndef MAKE_CONSENSUS
#define MAKE_CONSENSUS

#include <map>
#include <string>

class MakeConsensus {
public :
  MakeConsensus();
  
  void execute(std::map<std::string, std::string> & parameters);
};

#endif // MAKE_CONSENSUS
