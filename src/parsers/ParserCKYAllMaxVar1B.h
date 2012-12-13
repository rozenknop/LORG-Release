// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARONEBEST_H_
#define _PARSERCKYALLMAXVARONEBEST_H_

#include "edges/MaxRuleProbability1B.h"
#include "ParserCKYAllMaxVar.h"
#include "Word.h"


class ParserCKYAllMaxRule1B : public ParserCKYAllMaxRule<MaxRule1BTypes>
{
public:
  ParserCKYAllMaxRule1B(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn);

  ~ParserCKYAllMaxRule1B() {};

  void extract_solution();


protected:

  void change_rules_reset() const;


  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
     and uses this to select the best edge (max rule parsing)
   */
  void calculate_chart_specific_rule_probabilities_and_best_edge();

};

#endif /* _PARSERCKYALLMAXVARONEBEST_H_ */
