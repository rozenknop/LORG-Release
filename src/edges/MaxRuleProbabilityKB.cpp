// -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITYKB_CPP_
#define _MAXRULEPROBABILITYKB_CPP_

#include "MaxRuleProbabilityKB.h"


double MaxRuleProbabilityKB::log_normalisation_factor = 0;
unsigned MaxRuleProbabilityKB::size = 0;


std::ostream & MaxRuleProbabilityKB::operator>> (std::ostream & out) const
{
  for(auto& cand: candidates) { out << "cand:" << cand.probability << " "; }
  return out;
}




#endif /* _MAXRULEPROBABILITYKB_CPP_ */
