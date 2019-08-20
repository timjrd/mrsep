#pragma once

#include <cmath>
#include <tuple>
#include <vector>
#include <set>
#include <map>

#include "Base.hh"
#include "Instance.hh"

class Model {
private:
  std::vector<std::vector<Base>> _strainTypes;
  std::set<std::tuple<int,int,int>> _nullVars;
  std::map<std::tuple<int,int,int>,double> _weights;
  std::set<std::tuple<int,int,int>> _p;
  int _maxI = 0;
  std::map<int,int> _maxQ;
  int _maxK = 0;
  std::set<std::tuple<int,int>> _gap;
  std::map<int,int> _gaps;
  
  double const _epsilon;
  int const _maxMismatches;
  
  
  void addStrainType(std::vector<Base> const & strainType) {
    _strainTypes.push_back(strainType);
    int const j = _strainTypes.size();
    
    if (strainType.size() > _maxK)
      _maxK = strainType.size();

    int gaps = 0;
    for (int k = 1; k <= strainType.size(); k++)
      if (strainType[k-1] == Base::Empty) {
	gaps++;
	_gap.insert(std::make_tuple(j,k));
      }

    _gaps[j] = gaps;
  }

  // YOU MUST ADD ALL STRAIN TYPES BEFORE CALLING THIS METHOD.
  // i - read index
  // q - read mapping index
  // bases - nucleobases of the span with their associated positions
  void addReadMapping( int const i, int const q,
		       std::vector<std::tuple<Base,int>> const & bases ) {
    for (int j = 1; j <= maxJ(); j++) {
      int mismatches = 0;
      
      for (auto const & b : bases)
	if (_strainTypes[j-1].at(std::get<1>(b)-1) != std::get<0>(b))
	  mismatches++;

      if (mismatches <= _maxMismatches)
	// see Input::weight
	_weights[std::make_tuple(i,q,j)] = std::pow(0.5, mismatches);
      else
	// see Input::x
	_nullVars.insert(std::make_tuple(i,q,j));
    }
	
    for (auto const & b : bases) {
      int k = std::get<1>(b);
      if (k > _maxK)
	_maxK = k;
      // see Input::p
      _p.insert(std::make_tuple(i,q,k));
    }

    if (i > _maxI)
      _maxI = i;

    if (q > _maxQ[i])
      _maxQ[i] = q;
  }
  
  
public:
  // epsilon [0,1] - involved in constraints (3) and (4).
  Model( double const epsilon,
	 int const maxMismatches,
	 Instance const & instance )
    : _epsilon(epsilon)
    , _maxMismatches(maxMismatches) {
    
    // add strain types
    for (int i = 0; i < instance.positions().front.strainTypes.size(); i++) {
      std::vector<Base> strainType;
      
      for (auto const & p : instance.positions())
	strainType.push_back(p.strainTypes[i]);
      
      addStrainType(strainType);
    }
    
    // add read mappings
    for (int i = 1;; i++) {
      bool found = false;
      
      for (int q = 1;; q++) {
	std::vector<std::tuple<Base,int>> bases;
	
	int pos = 0;
	for (auto const & p : instance.positions()) {
	  for (auto const & rm : p.readMappings)
	    if (rm.read == i && rm.mapping == q)
	      bases.push_back({rm.base, pos});
	  
	  pos++;
	}
	
	if (not bases.empty()) {
	  addReadMapping(i, q, bases);
	  found = true;	  
	}
	else
	  break;
      }
      
      if (not found) break;
    }
  }
  
  double weight(std::tuple<int,int,int> const key) const {
    auto const it = _weights.find(key);
    if (it != _weights.end())
      return it->second;
    else
      return 0;
  }
  double weight(int const i, int const q, int const j) const {
    return weight(std::make_tuple(i,q,j));
  }
  
  // True if mapping q of read i is fully matching a subsequence of strain type j.
  bool x(int const i, int const q, int const j) const {
    return not _nullVars.count(std::make_tuple(i,q,j));
  }
  
  // True if mapping q of read i is covering position k / strain types column index k.
  bool p(int const i, int const q, int const k) const {
    return _p.count(std::make_tuple(i,q,k));
  }

  // Number of reads.
  int maxI() const {
    return _maxI;
  }

  // Number of mappings for read i.
  int maxQ(int const i) const {
    auto const it = _maxQ.find(i);
    if (it != _maxQ.end())
      return it->second;
    else
      return 0;
  }

  // Matrix width |L|.
  int maxK() const {
    return _maxK;
  }

  // Number of strain types.
  int maxJ() const {
    return _strainTypes.size();
  }

  // Is there a gap in strain type j at position k ?
  bool gap(int const j, int const k) const {
    return _gap.count(std::make_tuple(j,k));
  }

  // How many gaps are there in strain type j ?
  int gaps(int const j) const {
    auto const it = _gaps.find(j);
    if (it != _gaps.end())
      return it->second;
    else
      return 0;    
  }
  
  // [0,1] involved in constraints (3) and (4).
  double epsilon() const {
    return _epsilon;
  }
};
