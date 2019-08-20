#pragma once

#include <utility>
#include <tuple>
#include <set>
#include <map>
#include <ostream>

#include "Stats.hh"
#include "Instance.hh"

class Solver {
public:
  class Solution : public std::set<std::tuple<int,int,int>> {
  private:
    Instance const & _instance;
    
  public:
    Solution(Instance const & instance)
      : _instance(instance)
    {}
    
    void toStream(std::ostream & stream) const {
      // Each output line is formatted as "i q j" meaning mapping q of
      // read i is assigned to strain type j. The first line (Gurobi
      // license notification) should be ignored.
      for (auto const & assignment : *this)
	stream << std::get<0>(assignment) << " "
	       << std::get<1>(assignment) << " "
	       << _instance.positions().front().strainTypes.at(std::get<2>(assignment)-1).id+1 << "\n";
    }
  };
  
  virtual std::pair<Solution,Stats> solve(Instance const &) = 0;
};
