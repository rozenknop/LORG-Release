#include "parsers/ParserCKYAll.hpp"
#include "ParserCKYAllViterbi.h"
#include "grammars/GrammarAnnotated.hpp"

ParserCKYAllViterbi::ParserCKYAllViterbi(std::vector<AGrammar *>& cgs,
                                         const std::vector<double>& p, double b_t,
                                         const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants_,
                                         bool accurate_, unsigned min_beam, int stubborn)
: ParserCKYAll_Impl<ViterbiTypes>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn)
{
  fine_grammar = new ParserCKYAll::AGrammar(*cgs[cgs.size()-1]);
  fine_grammar->set_logmode(); // the Viterbi algorithms assume the fine grammar rules weights are log_probs

  //TODO maybe make this a parser option?
  //create the coarse-to-fine map

  std::vector<AGrammar*> all_grammars(cgs);
  all_grammars.push_back(fine_grammar);

  create_coarse_to_fine_mapping(all_grammars);

  Edge::set_unary_chains(fine_grammar->get_unary_decoding_paths());
}


void ParserCKYAllViterbi::extract_solution()
{
  const AnnotatedLabelsInfo& annotations_info = fine_grammar->get_annotations_info();


  chart->opencells_apply_bottom_up(
    [&annotations_info](Cell & cell){

      for(unsigned i = 0; i < cell.get_max_size(); ++i) {
        if(cell.exists_edge(i))  {
          cell.get_edge(i).get_prob_model().set_size(annotations_info.get_number_of_annotations(i));
          cell.get_edge(i).replace_rule_probabilities(0);
          cell.get_edge(i).apply(&ViterbiProbability::update_lexical);
          cell.get_edge(i).apply(&ViterbiProbability::update_binary);
        }
      }
    }
  );
}

inline
const ParserCKYAll::AGrammar& ParserCKYAllViterbi::get_fine_grammar() const
{
  return *fine_grammar;
}
