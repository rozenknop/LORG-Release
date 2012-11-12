// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVAR_H_
#define _PARSERCKYALLMAXVAR_H_

#include "ParserCKYAll.h"
#include "edges/PackedEdge.h"
#include "edges/maxrule_functions.h"
#include "edges/MaxRuleProbability.h"

class ParserCKYAllMaxRule;

// class MaxRulePackedEdgeData
// {
// private:
//   packed_edge_probability best;
// 
//   // fake methods for compatibility with maxrule_functions
// public:
//   const packed_edge_probability& get(unsigned) const
//   {return best;}
// 
//   bool has_solution(unsigned i) const {return i == 0;} ;
//   unsigned n_deriv() const {return 1;}
// 
//   friend class ParserCKYAllMaxRule;
// };



typedef PCKYAllCell< PackedEdge< MaxRuleProbability>> ParserCKYAllMaxRuleCell ;

class ParserCKYAllMaxRule : public ParserCKYAll_Impl<ParserCKYAllMaxRuleCell>
{
public:
  ParserCKYAllMaxRule(std::vector<AGrammar*>& cgs,
		      const std::vector<double>& p, double b_t,
		      const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants_,
		      bool accurate_, unsigned min_beam, int stubborn, unsigned cell_threads)
  : ParserCKYAll_Impl<ParserCKYAllMaxRuleCell>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn, cell_threads)
  {
    //TODO maybe make this a parser option?
    //create the coarse-to-fine map
    create_coarse_to_fine_mapping(grammars);

    Edge::set_viterbi_unary_chains(grammars[grammars.size() - 1]->get_unary_decoding_paths());
}

  ~ParserCKYAllMaxRule() {};

  void extract_solution();
protected:

  void change_rules_reset() const;


  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
     and uses this to select the best edge (max rule parsing)
   */
  void calculate_chart_specific_rule_probabilities_and_best_edge();
  
  static double log_normalisation_factor;
};

double ParserCKYAllMaxRule::log_normalisation_factor = 0;


  template<class Class, class Return, class... Args>
  auto myMethod( Class * object, Return(Class::*f)(Args... args)) -> function<Return(Args... args)> {
    return [object,f](Args... args){return (object->*f)(args...);};
  }

  
void ParserCKYAllMaxRule::extract_solution()
{
    
//   std::cout << *chart << std::endl;

  compute_inside_outside_probabilities();

//   std::cout << *chart << std::endl;

  calculate_chart_specific_rule_probabilities_and_best_edge();

//   std:e:cout << *chart << std::endl;

  // PCKYAllCell& root = chart->get_root();
  // if (!root.exists_edge(SymbolTable::instance_nt()->get_label_id(LorgConstants::tree_root_name)))
  //   std::cout << "no axiom at root" << std::endl;
}

void ParserCKYAllMaxRule::calculate_chart_specific_rule_probabilities_and_best_edge()
{
  double sentence_probability = std::log(get_sentence_probability());

  MaxRuleProbability::set_log_normalisation_factor(sentence_probability);

  chart->opencells_apply_bottom_up(
    [](Cell & cell){ cell.calculate_maxrule_probabilities(); }
  );
}

#endif /* _PARSERCKYALLMAXVAR_H_ */
