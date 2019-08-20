#pragma once

#include "Base.hh"
#include "Filter.hh"

class ConstFilter : public Filter {
public:
  Instance filter(Instance const & in) override {
    if ( not in.positions().empty()
	 and in.positions().front().strainTypes.size() == 1 ) {
      return in;
    }
    else {
      std::set<int> keep;
      
      for (int i = 0; i < in.positions().size(); i++) {
	auto const & p = in.positions()[i];	
	Base const first = p.strainTypes.front().base;
	
	for (auto const & st : p.strainTypes)
	  if (st.base != first) {
	    keep.insert(i);
	    break;
	  }
      }
      
      return in.filterPositions(keep);
    }
  }
};
