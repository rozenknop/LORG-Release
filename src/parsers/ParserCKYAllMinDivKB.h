// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARKB_H_
#define _PARSERCKYALLMAXVARKB_H_

#include "ParserCKYAll.h"
#include "edges/MaxRuleUpdater.h"

class MinDivProbabilityKB
{

public:

  typedef std::vector<packed_edge_probability_with_index> heap_type;
  typedef PackedEdge<MinDivProbabilityKB> Edge;
  typedef PCKYAllCell<Edge> Cell;
  typedef UnaryPackedEdgeDaughters<Cell> UnaryDaughters;
  typedef BinaryPackedEdgeDaughters<Cell> BinaryDaughters;
  typedef LexicalPackedEdgeDaughters LexicalDaughters;
  typedef MaxRuleUpdater<MinDivProbabilityKB> Updater;
  
private:

  heap_type candidates;
  heap_type derivations;


  static double log_normalisation_factor;
  static unsigned size;

public:

  MinDivProbabilityKB() :  candidates(), derivations() {candidates.reserve(50);};
  ~MinDivProbabilityKB() {};

  inline static void set_size(unsigned k) {size = k;}

  inline static void set_log_normalisation_factor(double lnf) {log_normalisation_factor = lnf;};

  inline const packed_edge_probability_with_index& get(unsigned idx) const {return derivations[idx];}
  inline packed_edge_probability& get(unsigned idx) { return derivations[idx]; }


  inline void update_lexical(Edge& e, const LexicalDaughters& dtr);
  inline void update_unary(Edge& e, const UnaryDaughters& dtr);
  inline void update_binary(Edge& e, const BinaryDaughters& dtr);
  inline void finalize();
  
  inline void find_succ(Edge*,packed_edge_probability_with_index& pep, bool licence_unaries);
  void extend_derivation(Edge*, unsigned, bool) ;

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





typedef PCKYAllCell<PackedEdge<MinDivProbabilityKB> > MinDivProbabilityCell ;


class ParserCKYAllMinDivKB : public ParserCKYAll_Impl<MinDivProbabilityCell>
{
 private:
  unsigned k;

 public:
  ParserCKYAllMinDivKB(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn, unsigned k_);

  ~ParserCKYAllMinDivKB() {};

  void extract_solution();

 private:
  void initialise_candidates();

  void extend_all_derivations();
};

#endif /* _PARSERCKYALLMAXVARKB_H_ */
