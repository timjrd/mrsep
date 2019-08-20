#pragma once

#include <map>

enum class Base { A, T, C, G, EMPTY };

std::map<char,Base> const CHAR_TO_BASE = {
  {'A', Base::A    },
  {'T', Base::T    },
  {'C', Base::C    },
  {'G', Base::G    },
  {'-', Base::EMPTY}
};

std::map<Base,char> const BASE_TO_CHAR = {
  {Base::A    , 'A'},
  {Base::T    , 'T'},
  {Base::C    , 'C'},
  {Base::G    , 'G'},
  {Base::EMPTY, '-'}
};
