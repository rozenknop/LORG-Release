// -*- mode: c++ -*-
#ifndef PARSERCKYALL_CPP
#define PARSERCKYALL_CPP

#include "ParserCKYAll.h"


template<>
GrammarAnnotated<BRuleC2f, URuleC2f, LexicalRuleC2f>::GrammarAnnotated(const std::string& filename)
    :
    Grammar<BRuleC2f, URuleC2f, LexicalRuleC2f>::Grammar(),
    AnnotatedContents(),
    viterbi_decoding_paths()
{
  std::map<short ,unsigned short> map;

  std::vector<BRule> bin;
  std::vector<URule> un;
  std::vector<LexicalRule> lex;


  BURuleInputParser::read_rulefile(filename, lex, un, bin, map, history_trees);

  label_annotations.set_num_annotations_map(map);

  lexical_rules.insert(lexical_rules.end(),lex.begin(),lex.end());
  unary_rules.insert(unary_rules.end(),un.begin(),un.end());
  binary_rules.insert(binary_rules.end(),bin.begin(),bin.end());


  // copied from void TrainingGrammar::uncompact_all_rules()
  for(std::vector<BRuleC2f>::iterator brule_it = binary_rules.begin();
      brule_it != binary_rules.end(); ++brule_it) {

    brule_it->uncompact(label_annotations.get_number_of_annotations(brule_it->get_rhs0()),
                        label_annotations.get_number_of_annotations(brule_it->get_rhs1()));
  }

  for(std::vector<URuleC2f>::iterator urule_it = unary_rules.begin();
      urule_it != unary_rules.end(); ++urule_it) {
    urule_it->uncompact(label_annotations.get_number_of_annotations(urule_it->get_rhs0()));
  }

  // lexical rules are not compacted
  // initialisation in  parserckyallfactory
}



ParserCKYAll::ParserCKYAll(std::vector<AGrammar*>& cgs,
                           const std::vector<double>& p,
                           double prior_threshold,
                           const annot_descendants_type& annot_descendants_,
                           bool accurate_,
                           unsigned min_beam, int stubborn)
    :
    Parser(cgs[0]),
    grammars(cgs),
    priors(p), prior_beam_threshold(prior_threshold),
    annot_descendants(annot_descendants_),
    accurate(accurate_),
    min_length_beam(min_beam),
    stubbornness(stubborn)
{
  // these thresholds look familiar ;)
  if(accurate) {
    //extra accurate
    //double t_acc[] = {-1000, -1000, -1000, -1000, -1000, -1000,-1000};
    double t_acc[] = {-8, -12, -12, -11, -12, -12, -14, -17};
    io_beam_thresholds = std::vector<double> (t_acc,t_acc+8);
  }
  else {
    double t[] ={-8, -9.75, -10, -9.6, -9.66, -8.01, -7.4, -10, -12};
    io_beam_thresholds = std::vector<double> (t,t+8);
  }
}



#endif /*PARSERCKYALL_CPP*/
