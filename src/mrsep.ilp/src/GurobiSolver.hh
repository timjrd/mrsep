#pragma once

#include <utility>
#include <tuple>
#include <set>
#include <map>

#include "gurobi_c++.h"

#include "Solver.hh"
#include "GurobiModel.hh"

class GurobiSolver : public Solver {
private:
  double const _epsilon;
  int    const _maxMismatches;
  
public:
  GurobiSolver(double const epsilon, int const maxMismatches)
    :_epsilon(epsilon), _maxMismatches(maxMismatches)
  {}
  
  std::pair<Solution,Stats> solve(Instance const & instance) override {
    GurobiModel const input(_epsilon, _maxMismatches, instance);
    
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);
    GRBModel model = GRBModel(env);

    /* Declaring decision variables. */
    std::map<std::tuple<int,int,int>,GRBVar> vars;
    // for each read i
    for (int i = 1; i <= input.maxI(); i++)
      // for each mapping q of read i
      for (int q = 1; q <= input.maxQ(i); q++)
	// for each strain type j
	for (int j = 1; j <= input.maxJ(); j++)
	  // if mapping q of read i is fully matching a subsequence of
	  // strain type j
	  if (input.x(i,q,j))
	    // declare a new binary decision variable (i,q,j) taking 1
	    // if mapping q of read i is assigned to strain type j, 0
	    // otherwise
	    vars[{i,q,j}] = model.addVar( 0, 1, 0, GRB_BINARY,
					  std::to_string(i) + "." +
					  std::to_string(q) + "." +
					  std::to_string(j) );
  
    int nbConstraints1NotNull = 0;
    /* Building constraints (1). */
    // for each read i
    for (int i = 1; i <= input.maxI(); i++) {
      bool null = true;
      GRBLinExpr sum = 0;

      // for each mapping q of read i
      for (int q = 1; q <= input.maxQ(i); q++)
	// for each strain type j
	for (int j = 1; j <= input.maxJ(); j++)
	  // when decision variable (i,q,j) exists / is not fixed to zero
	  for ( auto var = vars.find({i,q,j});
		var != vars.end();
		var  = vars.end() ) {
	    // add it to the sum
	    sum += var->second;
	    null = false;
	  }

      if (not null) {
	// declare a new constraint (1)
	model.addConstr(sum <= 1, "1." + std::to_string(i));
	nbConstraints1NotNull++;
      }
    }

    int nbCoveragesNotNull = 0;
    /* Building coverages. */
    std::map<std::tuple<int,int>,GRBLinExpr> cov;
    // for each strain type j
    for (int j = 1; j <= input.maxJ(); j++)
      // for each position k / strain types column index k
      for (int k = 1; k <= input.maxK(); k++)
	if (not input.gap(j,k)) {
	  bool null = true;
	  GRBLinExpr sum = 0;

	  // for each read i
	  for (int i = 1; i <= input.maxI(); i++)
	    // for each mapping q of read i
	    for (int q = 1; q <= input.maxQ(i); q++)
	      // if mapping q of read i is covering position k
	      if (input.p(i,q,k))
		// when decision variable (i,q,j) exists / is not fixed to zero
		for ( auto var = vars.find({i,q,j});
		      var != vars.end();
		      var  = vars.end() ) {
		  // add it to the sum
		  sum += var->second;
		  null = false;
		}
      
	  // coverage of strain type j at position k is the sum
	  cov[{j,k}] = sum;

	  if (not null)
	    nbCoveragesNotNull++;
	}
  
    int nbConstraints34 = 0;
    /* Building constraints (2), (3) and (4). */
    // for each strain type j
    for (int j = 1; j <= input.maxJ(); j++) {
      GRBLinExpr sum = 0;

      // for strain type j, add up coverages for each position k
      for (int k = 1; k <= input.maxK(); k++)
	if (not input.gap(j,k))
	  sum += cov.at({j,k});

      // make the average coverage of strain type j (2)
      GRBLinExpr const avg = sum / (input.maxK() - input.gaps(j));

      // for each position k
      for (int k = 1; k <= input.maxK(); k++)
	if (not input.gap(j,k)) {
	  std::string const suffix = std::to_string(j) + "." + std::to_string(k);
	  // declare a new constraint (3)
	  model.addConstr( cov.at({j,k}) >= (1-input.epsilon())*avg,
			   "3." + suffix );
	  // declare a new constraint (4)
	  model.addConstr( cov.at({j,k}) <= (1+input.epsilon())*avg,
			   "4." + suffix );
	  nbConstraints34++;
	}
    }
  
    /* Building objective function (6). */
    GRBLinExpr sum = 0;
    // add up all decision variables
    for (auto const & var : vars)
      sum += input.weight(var.first) * var.second;
    model.setObjective(sum, GRB_MAXIMIZE);
    
    //model.write("model-" + std::to_string(_epsilon) + ".lp");
    
    /* Optimizing. */
    model.optimize();
  
    /* Gathering the results. */
    Solution solution(instance);
    for (auto const & var : vars)
      if (var.second.get(GRB_DoubleAttr_X))
    	solution.insert(var.first);
    
    Stats stats;
    stats["nb_vars"]                   = vars.size();
    stats["nb_constraints_1_not_null"] = nbConstraints1NotNull;
    stats["nb_coverages_not_null"]     = nbCoveragesNotNull;
    stats["nb_constraints_3_4"]        = nbConstraints34;
    
    return {solution, stats};
  }
};
