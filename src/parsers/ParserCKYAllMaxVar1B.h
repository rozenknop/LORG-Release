// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_H_
#define _PARSERCKYALLMAXVARONEBEST_H_

#include "ParserCKYAllMaxVar.h"
#include "edges/MaxRuleProbability1B.h"


typedef PCKYAllCell< PackedEdge< MaxRuleProbability1B>> ParserCKYAllMaxRuleCell ;

class ParserCKYAllMaxRule1B : public ParserCKYAllMaxRule<ParserCKYAllMaxRuleCell>
{
public:
  ParserCKYAllMaxRule1B(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn)
      : ParserCKYAllMaxRule<ParserCKYAllMaxRuleCell>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn)
  {
    // this is not in the super class because maxn parsing uses a
    //different mapping
    //create the coarse-to-fine map
    this->create_coarse_to_fine_mapping(this->grammars);

    Edge::set_unary_chains(this->grammars[this->grammars.size() - 1]->get_unary_decoding_paths());
  }

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

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
