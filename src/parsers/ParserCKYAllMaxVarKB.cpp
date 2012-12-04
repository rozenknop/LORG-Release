// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARKB_CPP_
#define _PARSERCKYALLMAXVARKB_CPP_

#include "ParserCKYAllMaxVarKB.h"
#include "ParserCKYAll.hpp"
#include "utils/tick_count.h"

ParserCKYAllMaxRuleKB::ParserCKYAllMaxRuleKB(std::vector<AGrammar*>& cgs,
                                             const std::vector<double>& p, double b_t,
                                             const annot_descendants_type& annot_descendants_,
                                             bool accurate_, unsigned min_beam, int stubborn, unsigned k_)
: ParserCKYAllMaxRule<MaxRuleKBTypes>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn) , k(k_)
{
  // this is not in the super class because maxn parsing uses a
  //different mapping
  //create the coarse-to-fine map
  this->create_coarse_to_fine_mapping(this->grammars);

  Edge::set_unary_chains(this->grammars[this->grammars.size() - 1]->get_unary_decoding_paths());
}


void ParserCKYAllMaxRuleKB::extend_all_derivations()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  Cell& root = chart->get_root();

  if (!root.exists_edge(start_symbol))
    //   //   std::cout << "no axiom at root" << std::endl;
    return;

  for (unsigned i = 2; i <= k; ++i)
  {
    //      std::cout << "before extend" << std::endl;
    chart->get_root().get_edge(start_symbol).extend_derivation(i,true);
  }
}


void ParserCKYAllMaxRuleKB::extract_solution()
{
  {
    //     BLOCKTIMING("compute_inside_outside_probabilities");
    compute_inside_outside_probabilities();
  }
  {
//     BLOCKTIMING("initialize_candidates");
    initialise_candidates();
  }
  //  std::cout << "after init cand" << std::endl;

  {
//     BLOCKTIMING("extend_all_derivations");
    extend_all_derivations();
  }
}


void ParserCKYAllMaxRuleKB::initialise_candidates()
{

  double sentence_probability = std::log(get_sentence_probability());
  //  unsigned sent_size = chart->get_size();

  MaxRuleProbabilityKB::set_log_normalisation_factor(sentence_probability);
  MaxRuleProbabilityKB::set_size(k);

  calculate_maxrule_probabilities();
}





#endif /* _PARSERCKYALLMAXVARKB_CPP_ */
