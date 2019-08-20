#pragma once

#include "Instance.hh"

class Filter {
public:
  virtual Instance filter(Instance const &) = 0;
};
