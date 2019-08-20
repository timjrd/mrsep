#include <experimental/optional>
#include <iostream>
#include <fstream>

#include "Instance.hh"
#include "InstanceStats.hh"
#include "NullCovFilter.hh"
#include "ConstFilter.hh"
#include "Solver.hh"
#include "GurobiSolver.hh"
#include "GreedySolver.hh"

int main(int argc, char ** argv) {
  std::experimental::optional<Instance> instance = {};
  
  if (argc == 3) {
    if (std::string(argv[2]) == "-")
      instance = Instance::fromStream(std::cin);
    else
      instance = Instance::fromFile(argv[2]);
  }
  else {
    std::cerr << "USAGE: " << argv[0] << " epsilon strain_types_and_read_mappings.txt\n"
              << "       " << argv[0] << " epsilon -\n";
    return EXIT_FAILURE;
  }
  
  if (not instance) {
    std::cerr << "ERROR with input data.\n";
    return EXIT_FAILURE;    
  }
  
  double epsilon;
  try {
    epsilon = std::stod(std::string(argv[1]));
    if (epsilon < 0 or epsilon > 1)
      throw std::exception();
  }
  catch (std::exception const & e) {
    std::cerr << "ERROR with epsilon.\n";
    return EXIT_FAILURE;
  }
  
  InstanceStats(*instance).toStream("unfiltered.", std::cerr);
  Instance const filtered = ConstFilter().filter( NullCovFilter().filter(*instance) );
  InstanceStats(filtered).toStream("filtered.", std::cerr);
  
  GurobiSolver gurobi(epsilon, 2);
  GreedySolver greedy(0.975, 2);
  
  //Solver & solver = gurobi;
  Solver & solver = greedy;
  
  auto const solution = solver.solve(filtered);
  solution.second.toStream("solver.", std::cerr);
  solution.first.toStream(std::cout);
  
  return EXIT_SUCCESS;
}
