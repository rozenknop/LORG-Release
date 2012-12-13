// // -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITY1B_H_
#define _MAXRULEPROBABILITY1B_H_
#include <numeric>


#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"

class MaxRuleProbability1B;

struct MaxRule1BTypes {
  typedef MaxRuleProbability1B Best ;
  typedef emptystruct EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< MaxRule1BTypes > Edge ;
  typedef BasePackedEdge< MaxRule1BTypes > PEdge ;
  typedef UPackedEdge< MaxRule1BTypes > UEdge ;
  typedef LBPackedEdge< MaxRule1BTypes > LBEdge ;
  typedef PCKYAllCell< MaxRule1BTypes > Cell ;
  typedef ChartCKY< MaxRule1BTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MaxRule1BTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MaxRule1BTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MaxRule1BTypes> LexicalDaughter;
};


class MaxRuleProbability1B
{
public:
  typedef typename MaxRule1BTypes::Edge Edge;
  typedef typename MaxRule1BTypes::PEdge PEdge;
  typedef typename MaxRule1BTypes::UEdge UEdge;
  typedef typename MaxRule1BTypes::LBEdge LBEdge;
  typedef typename MaxRule1BTypes::Cell Cell;
  typedef typename MaxRule1BTypes::UnaryDaughter UnaryDaughter;
  typedef typename MaxRule1BTypes::BinaryDaughter BinaryDaughter;
  typedef typename MaxRule1BTypes::LexicalDaughter LexicalDaughter;
  typedef MaxRuleTreeLogProbaComputer<MaxRuleProbability1B> QInsideComputer;

private:
  packed_edge_probability best;
  static double log_normalisation_factor;
public:
  MaxRuleProbability1B() : best() {};
  ~MaxRuleProbability1B() {};

  static void set_log_normalisation_factor(double lnf);

  inline const packed_edge_probability& get(unsigned/*ignored*/) const {return best;}
  inline packed_edge_probability& get(unsigned /*ignored*/) {return best;}

  inline void update_lexical(LBEdge& e,  const LexicalDaughter& dtr);
  inline void update_unary(UEdge& e, const UnaryDaughter& dtr);
  inline void update_binary(LBEdge& e, const BinaryDaughter& dtr);
  inline void finalize();

  inline unsigned n_deriv() const {return 1;}

  inline bool has_solution(unsigned i) const {return i == 0;} ;
};

inline std::ostream& operator<<(std::ostream& out, const MaxRuleProbability1B & prob)
{
  return out << "((MaxRule1BProb: " << &prob << ")";
}





#endif /* _MAXRULEPROBABILITY_H_ */
