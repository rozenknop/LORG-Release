// -*- mode: c++ -*-

#include "ViterbiProbability.hpp"

std::ostream& operator<<(std::ostream& out, const ViterbiProbability & prob)
{
  return out << "((ViterbiProb: " << &prob << ")" << ")";
}


