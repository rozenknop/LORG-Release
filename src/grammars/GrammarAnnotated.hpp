// -*- mode: c++ -*-
#ifndef _GRAMMARANNOTATED_HPP_
#define _GRAMMARANNOTATED_HPP_

#include "GrammarAnnotated.h"
#include "grammars/Grammar.h"

template <typename Bin, typename Un, typename Lex>
GrammarAnnotated<Bin,Un,Lex>::GrammarAnnotated()
: 
Grammar<Bin,Un,Lex>(), 
AnnotatedContents(),
viterbi_decoding_paths()
{}

namespace {
  
  struct remove_unlikely_helper_struct
  {
    const double& threshold;
    remove_unlikely_helper_struct(const double& threshold_)
    : threshold(threshold_) {};
    
    template<class T>
    void operator() (T& rule) const {rule.remove_unlikely_annotations(threshold);}
  };


  
  ///  a structure that compares unary rules (according to their symbolic part)
  // A1 -> B1 (p1,...,pn) == A2 -> B2 (q1,...,qn) <=> A1 == A2 && B1 == B2
  // used when looking for unary chains
  template <typename T>
  struct compare_urule
  {
    const T& candidate;
    compare_urule(const T& r) : candidate(r) {}
    bool operator()(const T& test) const
    {
      return candidate.get_lhs() == test.get_lhs() && candidate.get_rhs0() == test.get_rhs0();}
  };



  template<typename Lex>
  std::vector<Lex> compute_new_unknowns(const std::vector<Lex>& lexical_rules)
  {
    typedef std::vector<Lex> LexVect;
    
    static int unknown_id = SymbolTable::instance_word().insert(LorgConstants::token_unknown);
    
    //map from tags to lexical rules
    uomap<int, Lex> tags2rules;
    
    std::set<int> tags_w_unk;
    
    for(typename LexVect::const_iterator i(lexical_rules.begin()); i != lexical_rules.end(); ++i) {
      if (i->get_rhs0() == unknown_id) {
        tags_w_unk.insert(i->get_lhs());
      }
    }
    
    
    for(typename LexVect::const_iterator i(lexical_rules.begin()); i != lexical_rules.end(); ++i) {
      
      const std::string& word = SymbolTable::instance_word().get_label_string(i->get_rhs0());
      
      
      //lexical to word signature
      if(!tags_w_unk.count(i->get_lhs()) && word.length() >= 3 && word.substr(0,3) == "UNK") {
        
        if(!tags2rules.count(i->get_lhs())) {
          Lex copy(*i);
          copy.set_rhs0(unknown_id);
          tags2rules.insert(std::make_pair(i->get_lhs(),copy));
        }
        else {
          
          Lex& old_rule = tags2rules.find(i->get_lhs())->second;
          std::vector<double>& old_probabilities = old_rule.get_probability();
          
          
          if(i->get_probability().size() > old_probabilities.size())
            old_probabilities.resize(i->get_probability().size());
          
          for(unsigned a = 0; a < i->get_probability().size(); ++a) {
            old_probabilities[a] += i->get_probability()[a];
          }
          
          
          
        }
      }
    }
    
    std::vector<Lex> result;
    for(typename uomap<int,Lex>::const_iterator i(tags2rules.begin()) ; i != tags2rules.end(); ++i) {
      result.push_back(i->second);
    }
    
    return result;
  }
}



template<typename Bin, typename Un, typename Lex>
void GrammarAnnotated<Bin,Un,Lex>::init()
{
  Grammar<Bin,Un,Lex>::init();

  // std::vector<Lex> nl = compute_new_unknowns(Grammar<Bin,Un,Lex>::lexical_rules);
  // Grammar<Bin,Un,Lex>::lexical_rules.insert(Grammar<Bin,Un,Lex>::lexical_rules.end(),nl.begin(),nl.end());

  //  std::clog << "uncompacting" << std::endl;

  for(typename std::vector<Un>::iterator urule_it(Grammar<Bin,Un,Lex>::unary_rules.begin());
      urule_it != Grammar<Bin,Un,Lex>::unary_rules.end(); ++urule_it) {
    urule_it->get_probability().resize(label_annotations.get_number_of_annotations(urule_it->get_lhs()));
    urule_it->uncompact(label_annotations.get_number_of_annotations(urule_it->get_rhs0()));
  }

  //  std::clog << "computing unary chains" << std::endl;

  compute_unary_chains(viterbi_decoding_paths,1, Grammar<Bin,Un,Lex>::unary_rules);

  // std::for_each(Grammar<Bin,Un,Lex>::unary_rules.begin(), Grammar<Bin,Un,Lex>::unary_rules.end(),
  //  		std::mem_fun_ref(&URule::compact));


  //  std::clog << "leaving init" << std::endl;

}

template<typename Bin, typename Un, typename Lex>
void GrammarAnnotated<Bin,Un,Lex>::compute_unary_chains(PathMatrix& decoding_paths,
                                                        unsigned max_path_length,
                                                        std::vector<Un>& unaries)
{
  bool update;

  unsigned path_length = 0;

  do {
    ++path_length;
    update = false;

    // vector of (Urule,lhs,matrix of annotations)
    // the rule is the newly created rule obtained by combining 2 unary rules
    // the lhs is the intermediate symbol (ie B in A->B + B->C)
    // the matrix: M_ij is the best annotation for B to go from A_i to C_j
    // if M_ij is -1 it means that the path ends here
    typedef std::pair<Un,std::pair<int,std::vector< std::vector<int> > > > rule2label_paths;
    std::vector<rule2label_paths> new_rules_with_info;

    for(typename std::vector<Un>::const_iterator it1(unaries.begin()); it1 != unaries.end(); ++it1) {
      const Un& rule1 = *it1;

      //      std::clog << "Rule1: " << &rule1 << "\t" << rule1 << std::endl;

      for(typename std::vector<Un>::const_iterator it2(unaries.begin()); it2 != unaries.end(); ++it2) {
        const Un& rule2 = *it2;

        //std::clog << "Rule2: " <<&rule2 << "\t" << rule2 << std::endl;

        //composing A->B with B->C to create A->C
        if(rule2.get_lhs() != rule1.get_rhs0())
          continue;

        //we don't want to create X->X
        if(rule2.get_rhs0() == rule1.get_lhs())
          continue;

        if(SymbolTable::instance_nt().get_label_string(rule1.get_rhs0())[0] == 't')
          continue;

        //        std::clog << "path_length: " << path_length << "\tpRule1: " << &rule1 << "\t" << rule1 << "\tRule2: " << &rule2 << "\t" << rule2 << std::endl;



        //composing probabilities
        std::vector< std::vector<double> > probabilities(label_annotations.get_number_of_annotations(rule1.get_lhs()),
                                                         std::vector<double>(label_annotations.get_number_of_annotations(rule2.get_rhs0()),0.0));
        //annotation path
        std::vector< std::vector<int> > apath(label_annotations.get_number_of_annotations(rule1.get_lhs()),
                                              std::vector<int>(label_annotations.get_number_of_annotations(rule2.get_rhs0()),-1));


        for(unsigned i = 0; i < probabilities.size(); ++i) {
          for(unsigned j = 0; j < probabilities[i].size(); ++j) {
            for(unsigned k = 0; k < rule1.get_probability()[i].size(); ++k) {

              double candidate = rule1.get_probability(i,k) * rule2.get_probability(k,j);
              // if(candidate > probabilities[i][j]) {
                //   probabilities[i][j] = candidate;
              //   apath[i][j] = k;
              // }
              probabilities[i][j] += candidate;
              apath[i][j]=0;
            }
          }
        }

        Un composition(URule(rule1.get_lhs(), rule2.get_rhs0(),probabilities));

        // std::clog << "From:\n\t" << rule1 << "\n and:\n\t" << rule2
        // 	  << "\nI created:\n\t" << composition << std::endl;


        new_rules_with_info.push_back(std::make_pair(composition,std::make_pair(rule1.get_rhs0(),apath)));

      }
    }


    //    std::clog << "updating probs" << std::endl;

    //updating probabilities
    for(typename std::vector<rule2label_paths>::const_iterator it(new_rules_with_info.begin()); it != new_rules_with_info.end(); ++it){

      const Un& rule_candidate = it->first;
      int intermediate_symb = it->second.first;
      const std::vector< std::vector<int> >& apath = it->second.second;

      // doess A->C already exists ?
      typename std::vector<Un>::iterator old_rule_iter = std::find_if(unaries.begin(),unaries.end(),compare_urule<Un>(rule_candidate));

      //no
      // we add the new rule and the new paths
      if(old_rule_iter == unaries.end()) {
        //        std::clog << "one new rule: " << rule_candidate << std::endl;
        unaries.push_back(rule_candidate);

        for(unsigned i = 0; i < apath.size(); ++i)
          for(unsigned j = 0; j < apath[i].size(); ++j)
            decoding_paths[std::make_pair(rule_candidate.get_lhs(),i)][std::make_pair(rule_candidate.get_rhs0(),j)] = std::make_pair(intermediate_symb,apath[i][j]);

          update = true;
      }
      //yes
      // we update probabilities (and paths)
      //if the new ones are better
      else {
        // for(unsigned i = 0; i < label_annotations.get_number_of_annotations(rule_candidate.get_lhs()); ++i)
        //   for(unsigned j = 0; j < label_annotations.get_number_of_annotations(rule_candidate.get_rhs0()); ++j) {
          //     if(old_rule_iter->get_probability(i,j) < rule_candidate.get_probability(i,j)) {
            //       old_rule_iter->set_probability(i,j,rule_candidate.get_probability(i,j));
            //       decoding_paths[std::make_pair(rule_candidate.get_lhs(),i)][std::make_pair(rule_candidate.get_rhs0(),j)] = std::make_pair(intermediate_symb,apath[i][j]);
            //       update = true;
            //     }
            //   }


            //        std::clog << "updating: " << rule_candidate << std::endl;

        double sum_old = 0.0;
        double sum_candidate = 0.0;
        for(unsigned i = 0; i < label_annotations.get_number_of_annotations(rule_candidate.get_lhs()); ++i)
          for(unsigned j = 0; j < label_annotations.get_number_of_annotations(rule_candidate.get_rhs0()); ++j) {
            sum_old += old_rule_iter->get_probability(i,j);
            sum_candidate += rule_candidate.get_probability(i,j);
          }

        if(sum_candidate > sum_old) {
          //	  std::clog << "old:\n" << *old_rule_iter << std::endl;
          old_rule_iter->get_probability() = rule_candidate.get_probability();
          //	  std::clog << "new:\n" << *old_rule_iter << std::endl;
          for(unsigned i = 0; i < label_annotations.get_number_of_annotations(rule_candidate.get_lhs()); ++i)
            for(unsigned j = 0; j < label_annotations.get_number_of_annotations(rule_candidate.get_rhs0()); ++j)
              decoding_paths[std::make_pair(rule_candidate.get_lhs(),i)][std::make_pair(rule_candidate.get_rhs0(),j)] = std::make_pair(intermediate_symb,apath[i][j]);
            update = true;
        }
      }
    }
  }
  while(update && path_length <= max_path_length);
}

template<typename Bin, typename Un, typename Lex>
void GrammarAnnotated<Bin,Un,Lex>::set_logmode()
{
  std::for_each(Grammar<Bin,Un,Lex>::binary_rules.begin(),
    Grammar<Bin,Un,Lex>::binary_rules.end(),
    std::mem_fun_ref(&Bin::set_logmode));
  std::for_each(Grammar<Bin,Un,Lex>::unary_rules.begin(),
    Grammar<Bin,Un,Lex>::unary_rules.end(),
    std::mem_fun_ref(&Un::set_logmode));
  std::for_each(Grammar<Bin,Un,Lex>::lexical_rules.begin(),
    Grammar<Bin,Un,Lex>::lexical_rules.end(),
    std::mem_fun_ref(&Lex::set_logmode));
}

template<typename Bin, typename Un, typename Lex>
void GrammarAnnotated<Bin,Un,Lex>::remove_unlikely_annotations_all_rules(double threshold)
{
  remove_unlikely_helper_struct remover(threshold);
  std::for_each(Grammar<Bin,Un,Lex>::binary_rules.begin(),
    Grammar<Bin,Un,Lex>::binary_rules.end(),
    remover);


  //  std::clog << "bin before: " << Grammar<Bin,Un,Lex>::binary_rules.size() << std::endl;
  Grammar<Bin,Un,Lex>::binary_rules.erase(std::remove_if(Grammar<Bin,Un,Lex>::binary_rules.begin(),
              Grammar<Bin,Un,Lex>::binary_rules.end(),
              std::mem_fun_ref(&Bin::is_empty)),
            Grammar<Bin,Un,Lex>::binary_rules.end());
  //  std::clog << "bin after: " << Grammar<Bin,Un,Lex>::binary_rules.size() << std::endl;


  std::for_each(Grammar<Bin,Un,Lex>::unary_rules.begin(),
    Grammar<Bin,Un,Lex>::unary_rules.end(),
    remover);

  //  std::clog << "un before: " << Grammar<Bin,Un,Lex>::unary_rules.size() << std::endl;
  Grammar<Bin,Un,Lex>::unary_rules.erase(std::remove_if(Grammar<Bin,Un,Lex>::unary_rules.begin(),
              Grammar<Bin,Un,Lex>::unary_rules.end(),
              std::mem_fun_ref(&Un::is_empty)),
          Grammar<Bin,Un,Lex>::unary_rules.end());
  //  std::clog << "un after: " << Grammar<Bin,Un,Lex>::unary_rules.size() << std::endl;


  std::for_each(Grammar<Bin,Un,Lex>::lexical_rules.begin(),
    Grammar<Bin,Un,Lex>::lexical_rules.end(),
    remover);

  //  std::clog << "lex before: " << Grammar<Bin,Un,Lex>::lexical_rules.size() << std::endl;
  Grammar<Bin,Un,Lex>::lexical_rules.erase(std::remove_if(Grammar<Bin,Un,Lex>::lexical_rules.begin(),
                Grammar<Bin,Un,Lex>::lexical_rules.end(),
                std::mem_fun_ref(&Lex::is_empty)),
            Grammar<Bin,Un,Lex>::lexical_rules.end());
  //  std::clog << "lex after: " << Grammar<Bin,Un,Lex>::lexical_rules.size() << std::endl;

}


template<typename Bin, typename Un, typename Lex>
GrammarAnnotated<Bin,Un,Lex> *
GrammarAnnotated<Bin,Un,Lex>::create_projection(const std::vector<std::vector<double> >& expected_counts,
            const std::vector<std::vector<std::vector<unsigned> > >& annotation_mapping) const
{

  std::vector<std::vector<double> > conditional_probabilities = calculate_conditional_probs(expected_counts, annotation_mapping);

  GrammarAnnotated * result =  new GrammarAnnotated(conditional_probabilities, annotation_mapping,
                                                    Grammar<Bin,Un,Lex>::binary_rules,
                                                    Grammar<Bin,Un,Lex>::unary_rules,
                                                    Grammar<Bin,Un,Lex>::lexical_rules);


  return result;
}


// compute probabilities p(\alpha X_i \beta | Y_j )
//in table  table[Y,j][X,i] or [Y,j][word_id,0]
template<typename Bin, typename Un, typename Lex>
void GrammarAnnotated<Bin,Un,Lex>::compute_transition_probabilities(uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >& result) const
{
  for (typename std::vector<Bin>::const_iterator b(Grammar<Bin,Un,Lex>::binary_rules.begin());
       b != Grammar<Bin,Un,Lex>::binary_rules.end(); ++b)
  {
    const std::vector<std::vector<std::vector<double> > >& probs = b->get_probability();

    uomap<unsigned, uomap<int, uomap<unsigned, double > > >& lhs_map = result[b->get_lhs()];


    for (unsigned i = 0; i < probs.size(); ++i)
    {

      uomap<unsigned, double >& rhs0_map = lhs_map[i][b->get_rhs0()];
      uomap<unsigned, double >& rhs1_map = lhs_map[i][b->get_rhs1()];
      const std::vector<std::vector<double> >& probs_i = probs[i];

      for (unsigned j = 0; j < probs_i.size(); ++j)
      {
        double& rhs0_entry = rhs0_map[j];
        const std::vector<double>& probs_i_j = probs_i[j];

        for (unsigned k = 0; k < probs_i_j.size(); ++k)
        {
          rhs0_entry += probs_i_j[k];
          rhs1_map[k] += probs_i_j[k];
        }
      }
    }
  }


  for (typename std::vector<Un>::const_iterator u(Grammar<Bin,Un,Lex>::unary_rules.begin());
       u != Grammar<Bin,Un,Lex>::unary_rules.end(); ++u)
  {
    const std::vector<std::vector<double> >& probs = u->get_probability();
    uomap<unsigned, uomap<int, uomap<unsigned, double > > >& lhs_map = result[u->get_lhs()];

    for (unsigned i = 0; i < probs.size(); ++i)
    {
      const std::vector<double>& probs_i = probs[i];
      uomap<unsigned,double>& rhs0_map = lhs_map[i][u->get_rhs0()];
      for (unsigned j = 0; j < probs_i.size(); ++j)
      {
        rhs0_map[j] += probs_i[j];
      }
    }
  }

  // for (typename std::vector<Lex>::const_iterator l(Grammar<Bin,Un,Lex>::lexical_rules.begin());
  //      l != Grammar<Bin,Un,Lex>::lexical_rules.end(); ++l)
  //   {
  //     const std::vector<double>& probs = l->get_probability();
  //     uomap<unsigned, uomap<int, uomap<unsigned, double > > >& lhs_map = result[l->get_lhs()];


  //     for (unsigned i = 0; i < probs.size(); ++i)
  //       {
  //         lhs_map[i][l->get_word()][0] += probs[i];
  //       }
  //   }
}

template<typename Bin, typename Un, typename Lex>
GrammarAnnotated<Bin,Un,Lex>::GrammarAnnotated(const std::vector<std::vector<double> >& conditional_probabilities,
                                               const std::vector<std::vector<std::vector<unsigned> > >& mapping,
                                               const std::vector<Bin>& old_binary_rules,
                                               const std::vector<Un>& old_unary_rules,
                                               const std::vector<Lex>& old_lexical_rules) :

                                               Grammar<Bin,Un,Lex>(), AnnotatedContents(), viterbi_decoding_paths()

{
  //set label annotations
  // TODO get rid of it ?
  std::vector<unsigned short> la(mapping.size());
  for (unsigned short i = 0; i < mapping.size(); ++i)
  {
    la[i] = mapping[i].size();
  }

  label_annotations.set_num_annotations_map(la);


  project_rule pr(conditional_probabilities, mapping);

  //  std::clog << "before binaries" << std::endl;
  Grammar<Bin,Un,Lex>::binary_rules.reserve(old_binary_rules.size());
  std::transform(old_binary_rules.begin(),old_binary_rules.end(),
                  std::back_inserter(Grammar<Bin,Un,Lex>::binary_rules),
                  pr);
  //  std::clog << "after binaries" << std::endl;
  // std::copy(Grammar<Bin,Un,Lex>::binary_rules.begin(),
  // 	    Grammar<Bin,Un,Lex>::binary_rules.end(),
  //   	    std::ostream_iterator<Bin> (std::clog, "\n"));

  //  std::clog << "before unaries" << std::endl;
  Grammar<Bin,Un,Lex>::unary_rules.reserve(old_unary_rules.size());
  std::transform(old_unary_rules.begin(),old_unary_rules.end(),
                  std::back_inserter(Grammar<Bin,Un,Lex>::unary_rules),
                  pr);
  //  std::clog << "after unaries" << std::endl;
  // std::copy(Grammar<Bin,Un,Lex>::unary_rules.begin(),
  // 	    Grammar<Bin,Un,Lex>::unary_rules.end(),
  //   	    std::ostream_iterator<Un> (std::clog, "\n"));

  //  std::clog << "before lexicals" << std::endl;
  Grammar<Bin,Un,Lex>::lexical_rules.reserve(old_lexical_rules.size());
  std::transform(old_lexical_rules.begin(),old_lexical_rules.end(),
                  std::back_inserter(Grammar<Bin,Un,Lex>::lexical_rules),
                  pr);
  //  std::clog << "after lexicals" << std::endl;

  // std::copy(Grammar<Bin,Un,Lex>::lexical_rules.begin(),
  // 	    Grammar<Bin,Un,Lex>::lexical_rules.end(),
  //   	    std::ostream_iterator<Lex> (std::clog, "\n"));

}


// assume only one annotation per symbol !!!!!
template<typename Bin, typename Un, typename Lex>
std::vector<double>
GrammarAnnotated<Bin,Un,Lex>::compute_priors() const
{
  uomap<int, uomap<unsigned, uomap<int, uomap< unsigned, double > > > > transition_probabilities;
  compute_transition_probabilities(transition_probabilities);

  std::vector<std::vector<double> > expected_counts;
  calculate_expected_counts(transition_probabilities, get_annotations_info(), expected_counts);

  // assume only one annotation !!!!!
  double sum = 0;
  for (unsigned i = 0; i < expected_counts.size(); ++i)
  {
    sum += expected_counts[i][0];
  }

  std::vector<double> res(expected_counts.size());

  for (unsigned i = 0; i < res.size(); ++i)
  {
    res[i] = expected_counts[i][0] / sum;
  }

  return res;
}



#endif /* _GRAMMARANNOTATED_HPP_ */
