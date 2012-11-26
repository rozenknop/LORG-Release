// -*- mode: c++ -*-
#ifndef _GRAMMARANNOTATED_CPP_
#define _GRAMMARANNOTATED_CPP_

#include "GrammarAnnotated.h"
#include <numeric>




// replace current mapping codomain with its image
void extend_mapping(std::vector< std::vector<unsigned> >& acc, const std::vector<std::vector<unsigned> >& nexts)
{
  if(acc.empty())
    acc = nexts;
  else {
    for(unsigned annot = 0; annot < acc.size(); ++annot) {
      std::vector<unsigned> new_results;
      for(std::vector<unsigned>::const_iterator i(acc[annot].begin()); i != acc[annot].end(); ++i)
        new_results.insert(new_results.end(), nexts[*i].begin(), nexts[*i].end());
      
      acc[annot] = new_results;
    }
  }
}


// compute annotation mapping between two levels
std::vector<std::vector<std::vector<unsigned> > > compute_mapping(unsigned from, unsigned to,
                                                                  const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants)
{
  std::vector<std::vector<std::vector<unsigned> > > result(SymbolTable::instance_nt().get_symbol_count());
  unsigned begin = from;
  unsigned end = to;
  
  while(begin < end) {
    
    const std::vector<std::vector< std::vector<unsigned> > >& mapping_begin = annot_descendants[begin];
    
    for(unsigned nt = 0; nt < mapping_begin.size(); ++nt)
      
      extend_mapping(result[nt], mapping_begin[nt]);
    
    ++begin;
  }
  
  return result;
}

std::vector<uomap<unsigned,unsigned> > invert_mapping(std::vector<std::vector<std::vector<unsigned> > > mapping)
{
  std::vector<uomap<unsigned,unsigned> > res(mapping.size());
  
  for (unsigned i = 0; i < res.size(); ++i)
  {
    for (unsigned j = 0; j < mapping[i].size(); ++j)
    {
      for (std::vector<unsigned>::const_iterator it(mapping[i][j].begin()); it != mapping[i][j].end(); ++it)
      {
        res[i][*it] = j;
      }
    }
  }
  
  return res;
}


std::vector<std::vector<double> >
calculate_conditional_probs(const std::vector<std::vector<double> >& expected_counts,
                            const std::vector<std::vector<std::vector<unsigned> > >& mapping)
{
  std::vector<std::vector<double> > result(expected_counts);
  
  for (unsigned i = 0; i < mapping.size(); ++i)
  {
    const std::vector<std::vector<unsigned> >& mapping_i = mapping[i];
    const std::vector<double>& exp_i = expected_counts[i];
    std::vector<double>& result_i = result[i];
    
    
    for (unsigned annot_i = 0; annot_i < mapping_i.size(); ++annot_i)
    {
      double sum = 0;
      const std::vector<unsigned>& m = mapping_i[annot_i];
      for (std::vector<unsigned>::const_iterator it(m.begin()); it != m.end(); ++it)
      {
        sum+= exp_i[*it];
      }
      for (std::vector<unsigned>::const_iterator it(m.begin()); it != m.end(); ++it)
      {
        result_i[*it] /= sum;
      }
    }
  }
  
  return result;
}
  


  LexicalRuleC2f project_rule::operator()(const LexicalRuleC2f& old_rule) const
  {
    int lhs = old_rule.get_lhs();
    LexicalRuleC2f new_rule(LexicalRule(lhs, old_rule.get_word(),
                                        std::vector<double>(mapping[lhs].size())));
    
    std::vector<double>& probs = new_rule.get_probability();
    const std::vector<double>& old_probs = old_rule.get_probability();
    
    const uomap<unsigned,unsigned>& inverted_lhs = inverted[lhs];
    const std::vector<double>& conditional_probabilities_lhs = conditional_probabilities[lhs];
    
    for (unsigned i = 0; i < old_probs.size(); ++i)
    {
      unsigned new_lhs_annot = inverted_lhs.find(i)->second;
      probs[new_lhs_annot] += old_probs[i] * conditional_probabilities_lhs[i];
    }
    
    return new_rule;
  }
  
  URuleC2f project_rule::operator()(const URuleC2f& old_rule) const
  {
    //    std::clog << old_rule << std::endl;
    
    int lhs = old_rule.get_lhs();
    int rhs = old_rule.get_rhs0();
    URuleC2f new_rule(URule(lhs, rhs,
                            std::vector< std::vector<double> >(mapping[lhs].size(),
                                                               std::vector<double>(mapping[rhs].size()))));
    
    const std::vector<std::vector<double> >& old_probs = old_rule.get_probability();
    std::vector<std::vector<double> >& new_probs = new_rule.get_probability();
    
    
    const uomap<unsigned,unsigned>& inverted_lhs = inverted[lhs];
    const uomap<unsigned,unsigned>& inverted_rhs = inverted[rhs];
    const std::vector<double>& conditional_probabilities_lhs = conditional_probabilities[lhs];
    
    
    
    for (unsigned i = 0; i < old_probs.size(); ++i)
    {
      
      const std::vector<double>& old_probs_i = old_probs[i];
      const double& conditional_probabilities_lhs_i = conditional_probabilities_lhs[i];
      
      std::vector<double>& new_probs_new_lhs_annot = new_probs[inverted_lhs.find(i)->second];
      
      for (unsigned j = 0; j < old_probs_i.size(); ++j)
      {
        unsigned new_rhs_annot = inverted_rhs.find(j)->second;
        new_probs_new_lhs_annot[new_rhs_annot] += old_probs_i[j] * conditional_probabilities_lhs_i;
      }
    }
    
    
    //    std::clog << new_rule << std::endl;
    
    return new_rule;
  }
  
  
  BRuleC2f project_rule::operator()(const BRuleC2f& old_rule) const
  {
    int lhs = old_rule.get_lhs();
    int rhs0 = old_rule.get_rhs0();
    int rhs1 = old_rule.get_rhs1();
    BRuleC2f new_rule(BRule(lhs, rhs0, rhs1,
                            std::vector< std::vector< std::vector<double> > >(mapping[lhs].size(),
                                                                              std::vector< std::vector<double> >(mapping[rhs0].size(),
                                                                                                                 std::vector<double>(mapping[rhs1].size())))));
    
    
    const std::vector< std::vector<std::vector<double> > >& old_probs = old_rule.get_probability();
    std::vector< std::vector<std::vector<double> > >& new_probs = new_rule.get_probability();
    
    const uomap<unsigned,unsigned>& lhs_map = inverted[lhs];
    const uomap<unsigned, unsigned>& rhs0_map = inverted[rhs0];
    const uomap<unsigned, unsigned>& rhs1_map = inverted[rhs1];
    
    for (unsigned i = 0; i < old_probs.size(); ++i)
    {
      std::vector< std::vector<double> >& new_probs_lhs = new_probs[lhs_map.find(i)->second];
      double conditional = conditional_probabilities[lhs][i];
      const std::vector< std::vector<double> >& old_probs_i = old_probs[i];
      
      for (unsigned j = 0; j < old_probs_i.size(); ++j)
      {
        std::vector<double>& new_probs_lhs_rhs0 = new_probs_lhs[rhs0_map.find(j)->second];
        const std::vector<double>& old_probs_i_j = old_probs_i[j];
        for (unsigned k = 0; k < old_probs_i_j.size(); ++k)
        {
          new_probs_lhs_rhs0[rhs1_map.find(k)->second] += old_probs_i_j[k] * conditional;
        }
      }
    }
    
    return new_rule;
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
GrammarAnnotated<Bin,Un,Lex>::GrammarAnnotated
(const std::vector<std::vector<double> >& conditional_probabilities,
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



#endif /* _GRAMMARANNOTATED_H_ */
