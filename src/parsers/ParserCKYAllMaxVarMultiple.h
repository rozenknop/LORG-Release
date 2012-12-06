// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARMULTIPLE_H_
#define _PARSERCKYALLMAXVARMULTIPLE_H_

#include "ParserCKYAllMaxVar.h"
#include "edges/MaxRuleProbabilityMultiple.h"


class ParserCKYAllMaxRuleMultiple : public ParserCKYAllMaxRule<MaxRuleMultipleTypes>
{
public:
  ParserCKYAllMaxRuleMultiple(std::vector<AGrammar*>& cgs,
                              const std::vector<double>& p, double b_t,
                              const std::vector< std::vector<AGrammar*> >& fgs,
                              const std::vector<annot_descendants_type>& all_annot_descendants_,
                              bool accurate_, unsigned min_beam, int stubborn, unsigned k);

  ~ParserCKYAllMaxRuleMultiple();

  /**
     \brief wraps the calculation of the best derivation
  */
  void extract_solution();


  const AGrammar& get_fine_grammar(unsigned i, unsigned j) const;

protected:
  /**
     \brief compute scores with all the fine grammars and back them up in the chart
  */
  void precompute_all_backups();

  void multiple_inside_outside_specific();


  /**
     \brief replace rules with their followers according to the defined mapping
     and reset annotations to zero (and resize thme to 1)
   */
  void change_rules_reset();

  /**
     \brief replace rules with their followers + size_grammar (to skip intermediate grammars)
     and replace current annotations with backed up ones at position backup_idx
  */
  void change_rules_load_backup(unsigned backup_idx, unsigned size_grammar) const;


  void modify_backup(unsigned backup_idx) const;

  /**
     \brief Calculates the chart specific rule probabilities of the packed edges in the chart
  */
  void calculate_maxrule_probabilities();

  /**
     \brief pick up the best derivation once the edge scores have been calculated
   */
  void calculate_best_edge();

  /**
     \brief for all edges in chart, backup current annotation
   */
  void backup_annotations() const;


protected: // attributes
  std::vector< std::vector<AGrammar*> >fine_grammars; ///< the additional grammars to be used to extract the solution
  std::vector<annot_descendants_type> all_annot_descendants; ///< all the annotations mapping for the grammars (base + fine ones)


private:
  unsigned nb_grammars;
  unsigned k;
  void initialise_candidates();
  void extend_all_derivations();
};

inline
const ParserCKYAll::AGrammar& ParserCKYAllMaxRuleMultiple::get_fine_grammar(unsigned i, unsigned j) const
{
  return *fine_grammars[i][j];
}



#endif /* _PARSERCKYALLMAXVARMULTIPLE_H_ */
