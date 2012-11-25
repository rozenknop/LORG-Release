// -*- mode: c++ -*-
#ifndef _PARSERCKYALLVITERBI_H_
#define _PARSERCKYALLVITERBI_H_

#include "ParserCKYAll.h"
#include "edges/ViterbiProbability.h"

typedef PCKYAllCell<PackedEdge<ViterbiProbability> > ParserCKYAllViterbiCell ;

class ParserCKYAllViterbi : public ParserCKYAll_Impl<ParserCKYAllViterbiCell>
{
public:
  ParserCKYAllViterbi(std::vector<AGrammar*>& cgs,
                      const std::vector<double>& p, double b_t,
                      const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants_, bool accurate_, unsigned min_beam, int stubborn);

  virtual ~ParserCKYAllViterbi() { delete fine_grammar; fine_grammar = NULL;};

  void extract_solution();
  const AGrammar& get_fine_grammar() const;

private: // attributes
  AGrammar* fine_grammar; ///< the grammar to be used to extract the solution
};

inline
const ParserCKYAll::AGrammar& ParserCKYAllViterbi::get_fine_grammar() const
{
  return *fine_grammar;
}


ParserCKYAllViterbi::ParserCKYAllViterbi(std::vector<AGrammar *>& cgs,
                                         const std::vector<double>& p, double b_t,
                                         const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants_,
                                         bool accurate_, unsigned min_beam, int stubborn)
: ParserCKYAll_Impl<ParserCKYAllViterbiCell>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn)
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



#endif /* _PARSERCKYALLVITERBI_H_ */
