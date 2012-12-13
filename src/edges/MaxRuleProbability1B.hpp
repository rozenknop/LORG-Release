// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY_HPP_
#define _MAXRULEPROBABILITY_HPP_

#include "MaxRuleProbability1B.h"
#include "edges/PackedEdge.hpp"

#include <numeric>


inline void MaxRuleProbability1B::update_lexical(LBEdge & edge,  const LexicalDaughter& dtr)
{
    const LexicalRuleC2f* rule = dtr.get_rule();
    
    //std::cout << "update with " << *rule << std::endl;
    
    double probability = QInsideComputer::compute(edge.get_annotations(), rule, log_normalisation_factor);
    
    if (probability > best.probability) {
        best.probability = probability;
        best.dtrs = &dtr;
    }
}


inline void MaxRuleProbability1B::update_unary (UEdge & e, const UnaryDaughter & dtr)
{
  double probability = -std::numeric_limits<double>::infinity();

  const LBEdge& left  = dtr.left_daughter();
  if(left.get_best().get(0).dtrs && (left.get_best().get(0).dtrs->is_lexical() || left.get_best().get(0).dtrs->is_binary())) {
    probability =  QInsideComputer::compute(e.get_annotations(), dtr, log_normalisation_factor);
  }

  if (probability > best.probability) {
    best.probability = probability;
    best.dtrs = &dtr;
  }
}

inline void
MaxRuleProbability1B::update_binary (LBEdge & e, const BinaryDaughter & dtr)
{
  double probability = QInsideComputer::compute(e.get_annotations(), dtr, log_normalisation_factor);

  if (probability > best.probability)
    {
      best.probability = probability;
      best.dtrs = &dtr;
    }
}

// //uncomment the function to get the best indexes
// // only useful if you want to print them out
inline void MaxRuleProbability1B::finalize() {}
// // {
// //   if(up) {
// //     //   std::cout << "here" << std::endl;
// //     if(best.dtrs !=  NULL) {
// //       if (best.dtrs->is_binary())
// //   {
// //     //    std::cout << "there" << std::endl;
// //     compute_best_indexes(up->get_annotations(),
// //              *static_cast<const BinaryPackedEdgeDaughters*>(best.dtrs),
// //              log_normalisation_factor,best.left_index,best.right_index);
// //   }
// //       else { //unary
// //   compute_best_indexes(up->get_annotations(),
// //            *static_cast<const UnaryPackedEdgeDaughters*>(best.dtrs),
// //            log_normalisation_factor,best.left_index);
// //       }
// //     }
// //   }
// // }







////////////


#endif /* _MAXRULEPROBABILITY_H_ */
