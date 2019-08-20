#pragma once

#include <string>
#include <map>
#include <ostream>

class Stats : public std::map<std::string,int> {
public:
  void toStream(std::string const & prefix, std::ostream & stream) const {
    for (auto const & stat : *this)
      stream << prefix << stat.first << " " << stat.second << "\n";    
  }
};
