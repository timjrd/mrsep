#pragma once

#include <algorithm>
#include <iostream>

#include <vector>
#include <map>

#include "Base.hh"
#include "Instance.hh"

class Coverage : public std::vector<std::pair<int,int>> {
public:
  static Coverage unambiguous(Instance const & inst, int const st) {
    Coverage result;
    
    for (int p = 0; p < inst.positions().size(); p++) {
      auto const & pos = inst.positions()[p];
      
      if (pos.strainTypes[st].base != Base::EMPTY) {
	int n = 0;
	
	for (auto const & r : pos.readMappings)
	  for (auto const & m : r.second)
	    if (inst.match(st, r.first, m.first)) {
	      bool ambiguous = false;
		
	      for (int st2 = 0; not ambiguous
		     and st2 < inst.positions().front().strainTypes.size();
		   st2++)
		if (st2 != st and inst.match(st2, r.first, m.first))
		  ambiguous = true;
		
	      if (not ambiguous)
		n++;
	    }

	result.push_back({n,p});
      }
    }

    std::sort(result.begin(), result.end());
    return result;
  }
  
  static Coverage selected( Instance const & inst, int const st,
			    std::map<int,std::set<int>> const & rms ) {
    Coverage result;
    
    for (int p = 0; p < inst.positions().size(); p++) {
      auto const & pos = inst.positions()[p];
      
      if (pos.strainTypes[st].base != Base::EMPTY) {
	int n = 0;
	
	for (auto const & r : pos.readMappings)
	  for (auto const & m : r.second) {
	    auto const q = rms.find(r.first);
	    if (q != rms.end() and q->second.count(m.first))
	      n++;
	  }
	
	result.push_back({n,p});
      }
    }
    
    std::sort(result.begin(), result.end());
    return result;    
  }
};

// for (int st = 0; st < _positions.front().strainTypes.size(); st++) {
//   int minPerfect     = INT_MAX;
//   int minUnambiguous = INT_MAX;
      
//   for (Position const & p : _positions)
//     if (p.strainTypes[st].base != Base::EMPTY) {
//       int perfect     = 0;
//       int unambiguous = 0;
	
//       for (auto const & r : p.readMappings)
// 	for (auto const & m : r.second)
// 	  if (match(st, r.first, m.first)) {
// 	    perfect++;

// 	    if (unambiguous < minUnambiguous) {
// 	      bool ambiguous = false;
// 	      for ( int st2 = 0;
// 		    not ambiguous and st2 < _positions.front().strainTypes.size();
// 		    st2++ )
// 		if (st2 != st and match(st2, r.first, m.first))
// 		  ambiguous = true;
		
// 	      if (not ambiguous)
// 		unambiguous++;
// 	    }
	      
// 	    if ( perfect >= minPerfect
// 		 and unambiguous >= minUnambiguous )
// 	      break;
// 	  }
	    
//       if (perfect < minPerfect)
// 	minPerfect = perfect;
	  
//       if (unambiguous < minUnambiguous)
// 	minUnambiguous = unambiguous;
//     }
      
//   std::string const id = std::to_string
//     ( _positions.front().strainTypes[st].id );
      
//   stats["min_perfect."     + id] = minPerfect;
//   stats["min_unambiguous." + id] = minUnambiguous;
// }
