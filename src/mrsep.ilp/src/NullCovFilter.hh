#pragma once

#include "Base.hh"
#include "Filter.hh"

class NullCovFilter : public Filter {
public:
  Instance filter(Instance const & in) override {
    std::set<int> keep;
    
    for (int st = 0; st < in.positions().front().strainTypes.size(); st++) {
      bool nullStrain = false;
      
      for (auto const & p : in.positions()) {
	if (p.strainTypes[st].base != Base::EMPTY) {
	  bool nullPos = true;
	  
	  for (auto const & r : p.readMappings)
	    for (auto const & m : r.second)
	      if (in.match(st, r.first, m.first)) {
		nullPos = false;
		break;
	      }
	  
	  if (nullPos) {
	    nullStrain = true;
	    break;
	  }
	}
      }
      
      if (not nullStrain)
	keep.insert(st);
    }
    
    return in.filterStrainTypes(keep);
  }
};
