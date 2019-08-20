#include <optional>
#include <iostream>
#include <fstream>

#include "Instance.hh"
#include "Model.hh"
#include "gurobiSolve.hh"

int main(int argc, char ** argv) {
  std::optional<Instance> instance = std::nullopt;
  
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

  if (instance == std::nullopt) {
    std::cerr << "ERROR in input data.\n";
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
  
  Model model(epsilon, 0, instance.value()); // -,2,-
  
  // try {
  // }
  // catch (std::out_of_range const & e) {
  // }
  
  Output results;
  std::map<std::string,int> stats;  
  gurobiSolve(model, results, stats);

  // Each output line is formatted as "i q j" meaning mapping q of
  // read i is assigned to strain type j. The first line (Gurobi
  // license notification) should be ignored.
  for (auto const & assignment : results)
    std::cout << std::get<0>(assignment) << " "
              << std::get<1>(assignment) << " "
              << std::get<2>(assignment) << "\n";

  for (auto const & stat : stats)
    std::cerr << stat.first << " " << stat.second << "\n";
  
  return EXIT_SUCCESS;
}
