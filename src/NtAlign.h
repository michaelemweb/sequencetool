#ifndef NT_ALIGN
#define NT_ALIGN

#include <map>
#include <string>

class NtAlign {
 public :
  NtAlign();
  
  void execute(std::map<std::string, std::string> & parameters);
};

#endif //NT_ALIGN
