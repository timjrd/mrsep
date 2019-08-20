#pragma once

#include <iostream>

#include <climits>
#include <utility>
#include <set>
#include <map>

#include "Solver.hh"

class GreedySolver : public Solver {
private:
  double const _minReadMappingWeight;
  int    const _coverageRank;
  
public:
  GreedySolver(double const minReadMappingWeight, int const coverageRank)
    : _minReadMappingWeight(minReadMappingWeight)
    , _coverageRank(coverageRank)
  {}
  
  std::map<std::tuple<int,int,int>, double>
  mkReadMappingWeights(Instance const & instance) const {
    
    std::map<std::tuple<int,int,int>, double> result;

    for (int st = 0; st < instance.positions().front().strainTypes.size(); st++)
      for (auto const & read : instance.readMappings())
	for (auto const & mapping : read.second) {
	  int matches = 0;
	    
	  for (auto const & pb : mapping.second)
	    if ( pb.second ==
		 instance.positions()[pb.first].strainTypes[st].base )
	      matches++;

	  double const w = double(matches) / double(mapping.second.size());

	  if (w >= _minReadMappingWeight)
	    result[{st, read.first, mapping.first}] = w;
	}

    return result;
  }
    
  std::pair<Solution,Stats> solve(Instance const & instance) override {

    // Gurobi is printing a mandatory line to stdout at startup, so we
    // do the same.
    std::cout << "greedy solver" << std::endl;
    
    Solution solution(instance);
    Stats    stats;
    
    std::set<int> strainTypes;
    std::set<int> reads;

    auto const readMappingWeights = mkReadMappingWeights(instance);
    
    for (int st = 0; st < instance.positions().front().strainTypes.size(); st++)
      strainTypes.insert(st);

    for (auto const & r : instance.readMappings())
      reads.insert(r.first);

    int prevReads = INT_MAX;
    
    while ( reads.size() < prevReads
	    and not reads.empty()
	    and not strainTypes.empty() ) {
      prevReads = reads.size();
      
      std::map<std::pair<int,int>,int> strainTypeBudgets;
      std::map<int,int>                strainTypeScores;
      std::map<int, std::map<int,int>> chosenMappings;
      
      for (int const strainType : strainTypes) {

	/* We choose the best mapping for each read. */

	std::map<int,std::set<int>> available;
	for (auto const & read : reads)
	  for (auto const & mapping : instance.readMappings().at(read)) {
	    auto const q = readMappingWeights.find({strainType, read, mapping.first});
	    if (q != readMappingWeights.end())
	      available[read].insert(mapping.first);
	  }

	int prevAvailable = INT_MAX;
	
	while ( available.size() < prevAvailable
		and not available.empty() ) {
	  prevAvailable = available.size();
	  
	  /* We compute the approximate coverage. */
	
	  Coverage const approximate =
	    Coverage::selected(instance, strainType, available);

	  /* We choose the best mapping for the next read. */
	
	  for (auto const & p : approximate) {
	    if (p.first > 0) {
	      for (auto const & read : instance.positions()[p.second].readMappings) {
		int    bestMapping = -1;
		double bestWeight  = -1;

		for (auto const & mapping : read.second) {
		  auto const q = available.find(read.first);
		  if (q != available.end() and q->second.count(mapping.first)) {

		    double const w = readMappingWeights
		      .at({strainType, read.first, mapping.first});

		    if (w > bestWeight) {
		      bestMapping = mapping.first;
		      bestWeight  = w;
		    }
		  }
		}

		if (bestMapping >= 0) {
		  chosenMappings[strainType][read.first] = bestMapping;
		  available.erase(read.first);
		  goto continue_;
		}
	      }
	    }
	  }
	  continue_:;
	}

	if (not chosenMappings.count(strainType))
	  strainTypeScores[strainType] = 0;
	
	else {
	  
	  /* We compute the coverage with chosen mappings only. */
	  
	  std::map<int,std::set<int>> chosen;
	  for (auto const & read : chosenMappings.at(strainType))
	    chosen[read.first].insert(read.second);

	  Coverage coverage = Coverage::selected(instance, strainType, chosen);

	  /* We compute the approximate quantity of this strain type. */

	  int const quantity = coverage[_coverageRank].first;
	
	  /* We compute the score of the strain type, and we cap the coverage. */

	  int score = 0;
	
	  for (auto & p : coverage) {
	    if (p.first > quantity)
	      p.first = quantity;
	  
	    score += p.first;
	    strainTypeBudgets[{strainType, p.second}] = p.first;
	  }
	
	  strainTypeScores[strainType] = score;
	}
      }

      /* We take the strain type with the highest score. */

      int bestStrainType = -1;
      int bestScore      = -1;
      
      for (auto const & s : strainTypeScores)
	if (s.second > bestScore) {
	  bestStrainType = s.first;
	  bestScore      = s.second;
	}

      /* We assign reads to the strain type. */

      auto const cm = chosenMappings.find(bestStrainType);
      if (cm != chosenMappings.end()) {
      
	std::map<int,std::set<int>> available;
	for (auto const & read : reads) {
	  auto const r = cm->second.find(read);
	  if (r != cm->second.end())
	    available[read].insert(r->second);
	}

	int prevReads2 = INT_MAX;
      
	while ( reads.size() < prevReads2
		and not reads.empty() ) {
	  prevReads2 = reads.size();

	  /* We recompute the coverage. */
	
	  Coverage coverage = Coverage::selected(instance, bestStrainType, available);

	  /* We chose the best next read. */

	  for (auto const & p : coverage) {
	    int & budget = strainTypeBudgets.at({bestStrainType, p.second});
	  
	    if (budget > 0) {
	      int    bestRead    = -1;
	      double bestWeight  = -1;

	      for (auto const & read : instance.positions()[p.second].readMappings) {
		if (available.count(read.first)) {
		
		  int const mapping = chosenMappings.at(bestStrainType).at(read.first);
		  double const w = readMappingWeights.at({bestStrainType, read.first, mapping});
		
		  if (w > bestWeight) {
		    bestRead    = read.first;
		    bestWeight  = w;
		  }
		}
	      }

	      if (bestRead >= 0) {

		budget--;
	      
		reads.erase(bestRead);
		available.erase(bestRead);

		solution.insert({ 1 + bestRead,
				  1 + chosenMappings.at(bestStrainType).at(bestRead),
				  1 + bestStrainType });
	      
		break;
	      }

	    
	    }
	  }
	}
      }

      strainTypes.erase(bestStrainType);
    }
    
    return {solution, stats};
  }
};
