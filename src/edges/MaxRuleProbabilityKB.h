// -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITYKB_H_
#define _MAXRULEPROBABILITYKB_H_

#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"


class MaxRuleProbabilityKB;

struct MaxRuleKBTypes {
  typedef MaxRuleProbabilityKB Best;
  typedef emptystruct EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< MaxRuleKBTypes > Edge ;
  typedef AnnotatedEdge< MaxRuleKBTypes > AEdge ;
  typedef BasePackedEdge< MaxRuleKBTypes > PEdge ;
  typedef UPackedEdge< MaxRuleKBTypes > UEdge ;
  typedef LBPackedEdge< MaxRuleKBTypes > LBEdge ;
  typedef PCKYAllCell< MaxRuleKBTypes > Cell ;
  typedef ChartCKY< MaxRuleKBTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MaxRuleKBTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MaxRuleKBTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MaxRuleKBTypes> LexicalDaughter;
};


class MaxRuleProbabilityKB
{

public:

  typedef std::vector<packed_edge_probability_with_index> heap_type;
  typedef typename MaxRuleKBTypes::Edge Edge;
  typedef typename MaxRuleKBTypes::AEdge AEdge;
  typedef typename MaxRuleKBTypes::PEdge PEdge;
  typedef typename MaxRuleKBTypes::UEdge UEdge;
  typedef typename MaxRuleKBTypes::LBEdge LBEdge;
  typedef typename MaxRuleKBTypes::Cell Cell;
  typedef typename MaxRuleKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MaxRuleKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MaxRuleKBTypes::LexicalDaughter LexicalDaughter;
  typedef MaxRuleTreeLogProbaComputer<MaxRuleProbabilityKB> QInsideComputer;
  
private:

  heap_type candidates;
  heap_type derivations;


  static double log_normalisation_factor;
  static unsigned size;

public:

  MaxRuleProbabilityKB() :  candidates(), derivations() {candidates.reserve(50);};
  ~MaxRuleProbabilityKB() {};

  inline static void set_size(unsigned k) {size = k;}

  inline static void set_log_normalisation_factor(double lnf) {log_normalisation_factor = lnf;};

  inline const heap_type & get_candidates() const { return candidates; }
  inline const heap_type & get_derivations() const { return derivations; }

  inline const packed_edge_probability_with_index& get(unsigned idx) const {return derivations[idx];}
  inline packed_edge_probability& get(unsigned idx) { return derivations[idx]; }


  inline void update_lexical(LBEdge& e, const LexicalDaughter& dtr);
  inline void update_unary(UEdge& e, const UnaryDaughter& dtr);
  inline void update_binary(LBEdge& e, const BinaryDaughter& dtr);
  inline void finalize();
  
  inline void find_succ(PEdge*,packed_edge_probability_with_index& pep, bool licence_unaries);
  inline void extend_derivation(PEdge*, unsigned, bool) ;

  inline unsigned n_deriv() const {return derivations.size();};

  inline bool has_solution(unsigned i) const {return i <derivations.size();}

private:
  
  struct test_helper
  {
    const packed_edge_probability_with_index& pep;
    test_helper(const packed_edge_probability_with_index& p) : pep(p) {};

    inline bool operator()(const packed_edge_probability_with_index& p)
    {
      return (p.probability == pep.probability) //|| (p.dtrs == pep.dtrs)
      ;
    }
  };
  
  public:
    inline std::ostream& operator>>(std::ostream& out) const;
};


inline std::ostream& operator<<(std::ostream& out, const MaxRuleProbabilityKB & prob)
{
  return out << "((MaxRuleKBProb: " << &prob << "): nb_deriv." << prob.get_derivations().size() << " nb_candid." << prob.get_candidates().size() << ")";
}



#endif /* _MAXRULEPROBABILITYKB_H_ */
