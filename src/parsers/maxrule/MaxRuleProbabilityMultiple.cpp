// -*- mode: c++ -*-

#ifndef _MAXRULEMULTIPLEPROBABILITY_CPP_
#define _MAXRULEMULTIPLEPROBABILITY_CPP_

#include "MaxRuleProbabilityMultiple.h"

unsigned MaxRuleProbabilityMultiple::size = 1;
double MaxRuleProbabilityMultiple::log_normalisation_factor = 0;
unsigned MaxRuleProbabilityMultiple::nb_grammars = 0;
std::vector<double> MaxRuleProbabilityMultiple::log_normalisation_factor_backup;


void MaxRuleProbabilityMultiple::set_log_normalisation_factor(double lnf)
{
  log_normalisation_factor = lnf;
  log_normalisation_factor_backup.push_back(lnf);
}

void MaxRuleProbabilityMultiple::reset_log_normalisation_factor()
{
  log_normalisation_factor_backup.resize(0);
}

const double& MaxRuleProbabilityMultiple::get_log_normalisation_factor()
{
  return log_normalisation_factor;
}

#endif /* _MAXRULEMULTIPLEPROBABILITY_CPP_ */
