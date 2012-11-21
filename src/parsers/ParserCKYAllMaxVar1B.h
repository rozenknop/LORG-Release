// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_H_
#define _PARSERCKYALLMAXVARONEBEST_H_

#include "ParserCKYAllMaxVar.h"
#include "edges/MaxRuleProbability.h"


typedef PCKYAllCell< PackedEdge< MaxRuleProbability>> ParserCKYAllMaxRuleCell ;

class ParserCKYAllMaxRule1B : public ParserCKYAllMaxRule<ParserCKYAllMaxRuleCell>
{
public:
  ParserCKYAllMaxRule1B(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn, unsigned cell_threads)
      : ParserCKYAllMaxRule<ParserCKYAllMaxRuleCell>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn, cell_threads)
  {}

  ~ParserCKYAllMaxRule1B() {};

  void extract_solution();


protected:

  void change_rules_reset() const;


  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
     and uses this to select the best edge (max rule parsing)
   */
  void calculate_chart_specific_rule_probabilities_and_best_edge();

  static double log_normalisation_factor;
};

#include "PCKYAllCell.h"
#include "edges/PackedEdge.h"
#include "edges/maxrule_functions.h"

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

  MaxRuleProbability::set_log_normalisation_factor(sentence_probability);
  ParserCKYAllMaxRule::calculate_maxrule_probabilities();
}

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
