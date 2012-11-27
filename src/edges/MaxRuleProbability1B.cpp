// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY_CPP_
#define _MAXRULEPROBABILITY_CPP_

#include "MaxRuleProbability1B.h"

#include <numeric>

double MaxRuleProbability::log_normalisation_factor = 0;

void MaxRuleProbability::set_log_normalisation_factor(double lnf)
{
  log_normalisation_factor = lnf;
}

#endif /* _MAXRULEPROBABILITY_H_ */
