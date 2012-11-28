// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMINDIVKB_H_
#define _PARSERCKYALLMINDIVKB_H_

#include "ParserCKYAll.h"
#include "edges/MaxRuleUpdater.h"
#include "emptystruct.h"

class MinDivProbabilityKB;

struct MinDivKBTypes {
  typedef MinDivProbabilityKB EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef PackedEdge< MinDivKBTypes > Edge ;
  typedef PCKYAllCell< MinDivKBTypes > Cell ;
  typedef ChartCKY< MinDivKBTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MinDivKBTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MinDivKBTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MinDivKBTypes> LexicalDaughter;
};

class MinDivProbabilityKB
{

public:

  typedef std::vector<packed_edge_probability_with_index> heap_type;

  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::Cell Cell;
  typedef typename MinDivKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MinDivKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MinDivKBTypes::LexicalDaughter LexicalDaughter;
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


  inline void update_lexical(Edge& e, const LexicalDaughter& dtr);
  inline void update_unary(Edge& e, const UnaryDaughter& dtr);
  inline void update_binary(Edge& e, const BinaryDaughter& dtr);
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





class ParserCKYAllMinDivKB : public ParserCKYAll_Impl<MinDivKBTypes>
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

#endif /* _PARSERCKYALLMINDIVKB_H_ */
