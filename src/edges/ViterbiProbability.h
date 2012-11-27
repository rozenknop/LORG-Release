// -*- mode: c++ -*-

#ifndef _VITERBIPROBABILITY_H_
#define _VITERBIPROBABILITY_H_

#include "PackedEdgeProbability.h"
#include "PackedEdge.h"


class ViterbiProbability //: public PackedEdgeProbability
{
public:
  typedef PackedEdge<ViterbiProbability> Edge;
  typedef PCKYAllCell<Edge> Cell;
  typedef UnaryPackedEdgeDaughters<Cell> UnaryDaughters;
  typedef BinaryPackedEdgeDaughters<Cell> BinaryDaughters;
  typedef LexicalPackedEdgeDaughters LexicalDaughters;
  
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

  void update_lexical(Edge& e,  const LexicalDaughters& dtr);
  void update_unary(Edge& e, const UnaryDaughters& dtr);
  void update_binary(Edge& e, const BinaryDaughters& dtr);
  void finalize() {};

  bool has_solution(unsigned i) const {return i == 0;} ;
};

#endif /* _VITERBIPROBABILITY_H_ */
