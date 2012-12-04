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

std::vector<uomap<unsigned,unsigned> > project_rule::invert_mapping(std::vector<std::vector<std::vector<unsigned> > > mapping)
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
  



 void calculate_expected_counts(uomap<int, uomap<unsigned, uomap<int, uomap<unsigned, double > > > >& trans,
                                const AnnotatedLabelsInfo& ali,
                                std::vector<std::vector<double> >& result)
{
  static int root_id = SymbolTable::instance_nt().get_label_id(LorgConstants::tree_root_name);
  unsigned n_nts = SymbolTable::instance_nt().get_symbol_count();

  result.resize(n_nts);
  result[root_id].resize(1);

  std::vector<std::vector<double> > temp(result);


  result[root_id][0] = 1;
  temp[root_id][0] = 1;

  unsigned n_iter = 30;
  double epsilon = 1;



  while(n_iter-- && epsilon > 1.0e-10) {

    //        std::clog << n_iter << " : " << epsilon << std::endl;

    for (size_t i = 0; i < result.size(); ++i)
      {
        // if( (i % 100) == 0) {
        //   std::clog << "i : " << i << "\r";
        //   std::flush(std::clog);
        // }

        std::vector<double>& result_i = result[i];
        //if(n_iter == 50)
          result_i.resize(ali.get_number_of_annotations(i));

        uomap< unsigned, uomap<int,uomap<unsigned,double> > >& map = trans[i];

        for (unsigned annot_i = 0; annot_i < result_i.size(); ++annot_i)
        {
          //std::clog << "annot_i : " << annot_i << std::endl;

          const double& old_value = result_i[annot_i];
          uomap<int,uomap<unsigned,double> >& lhs_map = map[annot_i];

          for (unsigned j = 0; j < result.size(); ++j)
          {

            //            std::clog << "j : " << j << std::endl;

            if((int) j == root_id) continue;

            std::vector<double>& temp_j = temp[j];
            //if(n_iter == 50)
            temp_j.resize(ali.get_number_of_annotations(j));

            uomap<int,uomap<unsigned,double> >::iterator found_key;
            if( (found_key = lhs_map.find(j)) != lhs_map.end() ) {

              //                std::clog << "here" << std::endl;
              for (unsigned annot_j = 0; annot_j < temp_j.size(); ++annot_j)
              {
                // std::clog << "annot_j : " << annot_j << std::endl;
                //                try {
                temp_j[annot_j] += old_value * found_key->second[annot_j];
                  //                }
                  //                catch(std::out_of_range& e) {}
                //                    std::clog << "there" << std::endl;
              }
            }
          }
        }
      }

    //    std::clog << "updated temp" << std::endl;

    epsilon = 0;

    for (unsigned i = 0; i < result.size(); ++i)
    {
      //        std::cout << SymbolTable::instance_nt().translate(i) << std::endl;

      if((int) i == root_id) {
        //          std::cout << result[0][0] << std::endl;
        continue;
      }

      std::vector<double>& temp_i = temp[i];
      std::vector<double>& result_i = result[i];

      for (unsigned annot_i = 0; annot_i < result_i.size(); ++annot_i)
      {
        epsilon += std::abs(temp_i[annot_i] - result_i[annot_i]);
        result_i[annot_i] = temp_i[annot_i];
        // std::cout << annot_i
        //           << " : "
        //           << result_i[annot_i]
        //           << std::endl;
        temp_i[annot_i] = 0;
      }
    }
  }
}
 

#endif /* _GRAMMARANNOTATED_H_ */
