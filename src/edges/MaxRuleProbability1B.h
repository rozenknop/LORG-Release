// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY1B_H_
#define _MAXRULEPROBABILITY1B_H_
#include <numeric>


#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "MaxRuleUpdater.h"


class MaxRuleProbability1B
{
public:
  typedef PackedEdge<MaxRuleProbability1B> Edge;
  typedef PCKYAllCell<Edge> Cell;
  typedef UnaryPackedEdgeDaughters<Cell> UnaryDaughters;
  typedef BinaryPackedEdgeDaughters<Cell> BinaryDaughters;
  typedef LexicalPackedEdgeDaughters LexicalDaughters;
  typedef MaxRuleUpdater<MaxRuleProbability1B> Updater;

private:
  packed_edge_probability best;
  static double log_normalisation_factor;
public:
  MaxRuleProbability1B() : best() {};
  ~MaxRuleProbability1B() {};

  static void set_log_normalisation_factor(double lnf);

  inline const packed_edge_probability& get(unsigned/*ignored*/) const {return best;}
  inline packed_edge_probability& get(unsigned /*ignored*/) {return best;}

  inline void update_lexical(Edge& e,  const LexicalDaughters& dtr);
  inline void update_unary(Edge& e, const UnaryDaughters& dtr);
  inline void update_binary(Edge& e, const BinaryDaughters& dtr);
  inline void finalize();

  inline unsigned n_deriv() const {return 1;}

  inline bool has_solution(unsigned i) const {return i == 0;} ;
};




inline void MaxRuleProbability1B::update_lexical(Edge & edge,  const LexicalDaughters& dtr)
{
    const LexicalRuleC2f* rule = dtr.get_rule();
    
    //std::cout << "update with " << *rule << std::endl;
    
    double probability = Updater::update_maxrule_probability(edge.get_annotations(), rule, log_normalisation_factor);
    
    if (probability > best.probability) {
        best.probability = probability;
        best.dtrs = &dtr;
    }
}


inline void MaxRuleProbability1B::update_unary (Edge & e, const UnaryDaughters & dtr)
{
  double probability = -std::numeric_limits<double>::infinity();

  Edge& left  = dtr.left_daughter()->get_edge(dtr.get_rule()->get_rhs0());
  if(left.get_prob_model().get(0).dtrs && (left.get_prob_model().get(0).dtrs->is_lexical() || left.get_prob_model().get(0).dtrs->is_binary())) {
    probability =  Updater::update_maxrule_probability(e.get_annotations(), dtr, log_normalisation_factor);
  }

  if (probability > best.probability) {
    best.probability = probability;
    best.dtrs = &dtr;
  }
}

inline void
MaxRuleProbability1B::update_binary (Edge & e, const BinaryDaughters & dtr)
{
  double probability = Updater::update_maxrule_probability(e.get_annotations(), dtr, log_normalisation_factor);

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
// // 	{
// // 	  //	  std::cout << "there" << std::endl;
// // 	  compute_best_indexes(up->get_annotations(),
// // 			       *static_cast<const BinaryPackedEdgeDaughters*>(best.dtrs),
// // 			       log_normalisation_factor,best.left_index,best.right_index);
// // 	}
// //       else { //unary
// // 	compute_best_indexes(up->get_annotations(),
// // 			     *static_cast<const UnaryPackedEdgeDaughters*>(best.dtrs),
// // 			     log_normalisation_factor,best.left_index);
// //       }
// //     }
// //   }
// // }







////////////



#endif /* _MAXRULEPROBABILITY_H_ */
