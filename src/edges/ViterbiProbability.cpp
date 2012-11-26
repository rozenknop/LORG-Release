// -*- mode: c++ -*-

#include "ViterbiProbability.h"

void ViterbiProbability::update_lexical(Edge& edge, const LexicalDaughters& dtr)
{
  const AnnotationInfo& a = edge.get_annotations();
  const LexicalRuleC2f* rule = dtr.get_rule();

  for (unsigned i = 0; i < rule->get_probability().size(); ++i) {
    if(a.valid_prob_at(i, LorgConstants::NullProba)) {
      double probability = rule->get_probability(i);
      packed_edge_probability& current_best = best[i];
      if (probability > current_best.probability) {
        best[i].probability = probability;
        current_best.dtrs = &dtr;
      }
    }
  }
}

void ViterbiProbability::update_unary(Edge& edge, const UnaryDaughters& dtr)
{
  const AnnotationInfo& a = edge.get_annotations();
  const std::vector<std::vector<double> >& rule_probs = dtr.get_rule()->get_probability();

  const PackedEdge<ViterbiProbability>& left = (dtr.left_daughter())->get_edge(dtr.get_rule()->get_rhs0());

  for (unsigned i = 0; i < rule_probs.size(); ++i) {
    if(!a.valid_prob_at(i, LorgConstants::NullProba)) continue;
    packed_edge_probability_with_index& current_best = best[i];
    const std::vector<double>& rule_probs_i = rule_probs[i];
    for (unsigned j = 0; j < rule_probs_i.size(); ++j) {
      if(!left.valid_prob_at(j)) continue;
      double probability = rule_probs_i[j] + left.get_prob_model().get(j).probability; // log-mode
      if (probability > current_best.probability) {
        current_best.probability = probability;
        current_best.dtrs = &dtr;
        current_best.left_index = j;
      }
    }
  }
}

void ViterbiProbability::update_binary(Edge& edge, const BinaryDaughters& dtr)
{
  const AnnotationInfo& a = edge.get_annotations();
  const std::vector<std::vector<std::vector<double> > >& rule_probs = dtr.get_rule()->get_probability();

  const PackedEdge<ViterbiProbability>& left  = dtr.left_daughter()->get_edge(dtr.get_rule()->get_rhs0());
  const PackedEdge<ViterbiProbability>& right = dtr.right_daughter()->get_edge(dtr.get_rule()->get_rhs1());

  for (unsigned i = 0; i < rule_probs.size(); ++i) {
    if(!a.valid_prob_at(i, LorgConstants::NullProba)) continue;
    packed_edge_probability_with_index& current_best = best[i];
    const std::vector<std::vector<double> >& rule_probs_i = rule_probs[i];
    for (unsigned j=0; j < rule_probs_i.size(); ++j) {
      if(!left.valid_prob_at(j)) continue;
      const std::vector<double>& rule_probs_ij = rule_probs_i[j];
      const double& left_best = left.get_prob_model().get(j).probability;
      for (unsigned k = 0; k < rule_probs_ij.size(); ++k) {
        if(!right.valid_prob_at(k)) continue;
        double probability = rule_probs_ij[k] + left_best + right.get_prob_model().get(k).probability; //log
        if (probability > current_best.probability) {
          //std::cout << " best so far " << std::endl;
          current_best.probability = probability;
          current_best.dtrs = &dtr;
          current_best.left_index = j;
          current_best.right_index= k;
        }
      }
    }
  }
}



#endif /* _VITERBIPROBABILITY_H_ */
