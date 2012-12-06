// -*- mode: c++ -*-

#ifndef _VITERBIPROBABILITY_H_
#define _VITERBIPROBABILITY_H_

#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "ChartCKY.h"
#include "emptystruct.h"

class ViterbiProbability;

struct ViterbiTypes {
  typedef ViterbiProbability EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< ViterbiTypes > Edge ;
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
  typedef typename ViterbiTypes::Cell Cell;
  typedef typename ViterbiTypes::UnaryDaughter  UnaryDaughter;
  typedef typename ViterbiTypes::BinaryDaughter BinaryDaughter;
  typedef typename ViterbiTypes::LexicalDaughter LexicalDaughter;
  
private:
  std::vector<packed_edge_probability_with_index> best;
public:
  ViterbiProbability() {};
  ViterbiProbability(unsigned size) : best(size) {};

  void set_size(unsigned size_) {best.resize(size_);}


  const packed_edge_probability& get(unsigned index) const
  {return best[index];}

  packed_edge_probability& get(unsigned index)
  {return best[index];}

  void update_lexical(Edge& e,  const LexicalDaughter& dtr);
  void update_unary(Edge& e, const UnaryDaughter& dtr);
  void update_binary(Edge& e, const BinaryDaughter& dtr);
  void finalize() {};

  bool has_solution(unsigned i) const {return i == 0;} ;
};

std::ostream& operator<<(std::ostream& out, const ViterbiProbability & prob);

#endif /* _VITERBIPROBABILITY_H_ */
