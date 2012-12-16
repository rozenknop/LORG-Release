// -*- mode: c++ -*-
#ifndef VITERBI_PROBABILITY_HPP
#define VITERBI_PROBABILITY_HPP

#include "ViterbiProbability.h"
#include "edges/PackedEdge.hpp"


inline void ViterbiProbability::update_lexical(LBEdge& edge, const LexicalDaughter& dtr)
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

inline void ViterbiProbability::update_unary(UEdge& edge, const UnaryDaughter& dtr)
{
  const AnnotationInfo& a = edge.get_annotations();
  const std::vector<std::vector<double> >& rule_probs = dtr.get_rule()->get_probability();

  const LBEdge& left = dtr.left_daughter();

  for (unsigned i = 0; i < rule_probs.size(); ++i) {
    if(!a.valid_prob_at(i, LorgConstants::NullProba)) continue;
    packed_edge_probability_with_index& current_best = best[i];
    const std::vector<double>& rule_probs_i = rule_probs[i];
    for (unsigned j = 0; j < rule_probs_i.size(); ++j) {
      if(!left.valid_prob_at(j)) continue;
      double probability = rule_probs_i[j] + left.get_best().get(j).probability; // log-mode
      if (probability > current_best.probability) {
        current_best.probability = probability;
        current_best.dtrs = &dtr;
        current_best.left_index = j;
      }
    }
  }
}

inline void ViterbiProbability::update_binary(LBEdge& edge, const BinaryDaughter& dtr)
{
  const AnnotationInfo& a = edge.get_annotations();
  const std::vector<std::vector<std::vector<double> > >& rule_probs = dtr.get_rule()->get_probability();

  const PEdge& left  = (PEdge &) dtr.left_daughter();
  const PEdge& right = (PEdge &) dtr.right_daughter();

  for (unsigned i = 0; i < rule_probs.size(); ++i) {
    if(!a.valid_prob_at(i, LorgConstants::NullProba)) continue;
    packed_edge_probability_with_index& current_best = best[i];
    const std::vector<std::vector<double> >& rule_probs_i = rule_probs[i];
    for (unsigned j=0; j < rule_probs_i.size(); ++j) {
      if(!left.valid_prob_at(j)) continue;
      const std::vector<double>& rule_probs_ij = rule_probs_i[j];
      const double& left_best = left.get_best().get(j).probability;
      for (unsigned k = 0; k < rule_probs_ij.size(); ++k) {
        if(!right.valid_prob_at(k)) continue;
        double probability = rule_probs_ij[k] + left_best + right.get_best().get(k).probability; //log
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

#endif
