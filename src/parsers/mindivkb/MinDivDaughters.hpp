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
  mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                              leftannot.inside_probabilities.array,
                                                              rightannot.inside_probabilities.array,
                                                              leftannot.outside_probabilities.array,
                                                              rightannot.outside_probabilities.array)
  / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
}

inline double MinDivBinaryDaughter::tree_log_proba(unsigned left_idx, unsigned right_idx) const
{
  return 
  log(q) 
  + left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx).probability
  + right_daughter()->get_edge(get_rule()->get_rhs1()).get_prob_model().get(right_idx).probability;
}

/***********************************************************/
/*              Unary daughters                            */
/***********************************************************/
inline void MinDivUnaryDaughter::outside_and_marginal(AnnotationInfo & annotations)
{
  auto & leftannot = left_daughter()->get_edge(get_rule()->get_rhs0()).get_annotations();
  mp = RH::rule->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                            leftannot.inside_probabilities.array,
                                                            leftannot.outside_probabilities_unary_temp.array)
        / MinDivProbabilityKB::get_normalisation_factor();
}
inline double MinDivUnaryDaughter::tree_log_proba(unsigned left_idx) const {
  return 
    log(q) 
    + left_daughter()->get_edge(get_rule()->get_rhs0()).get_prob_model().get(left_idx).probability;
}


  
/***********************************************************/
/*              Lexical daughters                          */
/***********************************************************/
inline void MinDivLexicalDaughter::outside_and_marginal(AnnotationInfo & annotations)
{
  mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array)
  / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
}

inline double MinDivLexicalDaughter::tree_log_proba() const { return log(q); }


#endif
