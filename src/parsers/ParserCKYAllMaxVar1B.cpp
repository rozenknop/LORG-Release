// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_CPP_
#define _PARSERCKYALLMAXVARONEBEST_CPP_

#include "ParserCKYAllMaxVar1B.h"
#include "ParserCKYAll.hpp"

double ParserCKYAllMaxRule1B::log_normalisation_factor = 0;

void ParserCKYAllMaxRule1B::extract_solution()
{
  //   std::cout << *chart << std::endl;
  compute_inside_outside_probabilities();
  //   std::cout << *chart << std::endl;
  calculate_chart_specific_rule_probabilities_and_best_edge();
  //   std:e:cout << *chart << std::endl;
  // PCKYAllCell& root = chart->get_root();
  // if (!root.exists_edge(SymbolTable::instance_nt()->get_label_id(LorgConstants::tree_root_name)))
  //   std::cout << "no axiom at root" << std::endl;
}

void ParserCKYAllMaxRule1B::calculate_chart_specific_rule_probabilities_and_best_edge()
{
  double sentence_probability = std::log(get_sentence_probability());
  MaxRuleProbability1B::set_log_normalisation_factor(sentence_probability);
  ParserCKYAllMaxRule::calculate_maxrule_probabilities();
}

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
