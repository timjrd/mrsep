#pragma once

#include "Stats.hh"
#include "Instance.hh"
#include "Coverage.hh"

class InstanceStats : public Stats {
public:
  InstanceStats(Instance const & inst) {
    (*this)["nb_positions"]    = inst.positions().size();
    (*this)["nb_strain_types"] = inst.positions().front().strainTypes.size();
    
    int n = 0;
    for (auto const & r : inst.readMappings())
      n += r.second.size();
    (*this)["nb_read_mappings"] = n;

    // for (int st = 0; st < inst.positions().front().strainTypes.size(); st++) {
    //   int const id = inst.positions().front().strainTypes[st].id;
    //   Coverage const cov = Coverage::unambiguous(inst, st);
    //   //for (int k = 0; k <= 150; k += 10)
    //   (*this)[ "min_unambiguous." +
    // 	       //std::to_string(k) + "." +
    // 	       std::to_string(id) ] = cov[150].first;
    // }
  }
};
