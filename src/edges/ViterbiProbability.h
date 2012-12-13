// -*- mode: c++ -*-

#ifndef _VITERBIPROBABILITY_H_
#define _VITERBIPROBABILITY_H_

#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "ChartCKY.h"
#include "emptystruct.h"

class ViterbiProbability;

struct ViterbiTypes {
  typedef ViterbiProbability Best ;
  typedef emptystruct EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< ViterbiTypes > Edge ;
  typedef BasePackedEdge<ViterbiTypes> PEdge;
  typedef UPackedEdge<ViterbiTypes> UEdge;
  typedef LBPackedEdge<ViterbiTypes> LBEdge;
  typedef PCKYAllCell< ViterbiTypes > Cell ;
  typedef ChartCKY< ViterbiTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<ViterbiTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<ViterbiTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<ViterbiTypes> LexicalDaughter;
};


class ViterbiProbability //: public PackedEdgeProbability
{
public:
  typedef typename ViterbiTypes::Edge Edge;
  typedef typename ViterbiTypes::PEdge PEdge;
  typedef UPackedEdge<ViterbiTypes> UEdge;
  typedef LBPackedEdge<ViterbiTypes> LBEdge;
  typedef typename ViterbiTypes::Cell Cell;
  typedef typename ViterbiTypes::UnaryDaughter  UnaryDaughter;
  typedef typename ViterbiTypes::BinaryDaughter BinaryDaughter;
  typedef typename ViterbiTypes::LexicalDaughter LexicalDaughter;
  
private:
  std::vector<packed_edge_probability_with_index> best;
public:
  ViterbiProbability() {};
  ViterbiProbability(unsigned size) : best(size) {};

  inline void set_size(unsigned size_) {best.resize(size_);}


  inline const packed_edge_probability& get(unsigned index) const
  {return best[index];}

  inline packed_edge_probability& get(unsigned index)
  {return best[index];}

  inline void update_lexical(LBEdge& e,  const LexicalDaughter& dtr);
  inline void update_unary(UEdge& e, const UnaryDaughter& dtr);
  inline void update_binary(LBEdge& e, const BinaryDaughter& dtr);
  inline void finalize() {};

  inline bool has_solution(unsigned i) const {return i == 0;} ;
};

std::ostream& operator<<(std::ostream& out, const ViterbiProbability & prob);

#endif /* _VITERBIPROBABILITY_H_ */
