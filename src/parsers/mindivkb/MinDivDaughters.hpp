#ifndef _MINDIVDAUGHTERS_HPP_
#define _MINDIVDAUGHTERS_HPP_

#include "MinDivDaughters.h"



/***********************************************************/
/*              Binary daughters                           */
/***********************************************************/
inline void MinDivBinaryDaughter::outside_and_marginal(AnnotationInfo & annotations)
{
  auto & leftannot = left_daughter()->get_edge(get_rule()->get_rhs0()).get_annotations();    
  auto & rightannot= right_daughter()->get_edge(get_rule()->get_rhs1()).get_annotations();
  q = mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                              leftannot.inside_probabilities.array,
                                                              rightannot.inside_probabilities.array,
                                                              leftannot.outside_probabilities.array,
                                                              rightannot.outside_probabilities.array)
  / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
}

inline double MinDivBinaryDaughter::tree_log_proba(unsigned left_idx, unsigned right_idx) const
{
  return std::min(0.0,
                  log(q) 
                  + left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx).probability
                  + right_daughter()->get_edge(get_rule()->get_rhs1()).get_prob_model().get(right_idx).probability);
}

/***********************************************************/
/*              Unary daughters                            */
/***********************************************************/
inline void MinDivUnaryDaughter::outside_and_marginal(AnnotationInfo & annotations)
{
  auto & leftannot = left_daughter()->get_edge(get_rule()->get_rhs0()).get_annotations();
  q = mp = RH::rule->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                                leftannot.inside_probabilities.array,
                                                                leftannot.outside_probabilities_unary_temp.array)
  / MinDivProbabilityKB::get_normalisation_factor();
}
inline double MinDivUnaryDaughter::tree_log_proba(unsigned left_idx) const {
  #warning dirty hack, false result
  // the chart should not aggregate heads of unary rules with heads of binary and lexical rules
  // see edges/MaxRuleTreeLogProbaComputer.h:195
  if (left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().n_deriv() > left_idx
    /*and left_daughter()->get_edge(get_rule()->get_rhs0()).get_unary_daughters().size() == 0*/)
  {
//     std::cout << "left " << left << ": " << *left << std::endl;
//     std::cout << "a_" << get_rule() << std::endl;
//     std::cout << "b_" << *get_rule() << std::endl;
//     std::cout << "c_" << get_rule()->get_rhs0() << std::endl;
//     std::cout << "d_" <<  &left_daughter()->get_edge(get_rule()->get_rhs0()) << std::endl;
//     std::cout << "e_" << & left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model() << std::endl;
//     std::cout << "f_" <<  &left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx) << std::endl;
//     std::cout << "g_" <<  left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx).probability << std::endl;
    
    return 
    std::min(0.0, 
             log(q) 
             + left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx).probability);
  }
  return -std::numeric_limits < double >::infinity ();
  }


  
/***********************************************************/
/*              Lexical daughters                          */
/***********************************************************/
inline void MinDivLexicalDaughter::outside_and_marginal(AnnotationInfo & annotations)
{
  q = mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array)
  / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
}

inline double MinDivLexicalDaughter::tree_log_proba() const {
  return std::min(0.0, log(q));
}


#endif
