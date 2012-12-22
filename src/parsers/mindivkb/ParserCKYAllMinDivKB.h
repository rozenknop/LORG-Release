// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMINDIVKB_H_
#define _PARSERCKYALLMINDIVKB_H_

#include "parsers/ParserCKYAll.h"
#include "emptystruct.h"
#include "MinDivTypes.h"


class MinDivBest
{
public:
  MinDivBest() { candidates.reserve(50); }

private:
  friend class ParserCKYAllMinDivKB;
  
  typedef std::vector<packed_edge_probability_with_index> heap_type;
  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::PEdge PEdge;
  typedef typename MinDivKBTypes::UEdge UEdge;
  typedef typename MinDivKBTypes::LBEdge LBEdge;
  typedef typename MinDivKBTypes::Cell Cell;
  typedef typename MinDivKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MinDivKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MinDivKBTypes::LexicalDaughter LexicalDaughter;

  heap_type candidates;
  heap_type derivations;
  static unsigned size;
  
public:
  inline static void set_size(unsigned k) {size = k;}
  
  inline const heap_type & get_candidates() const { return candidates; }
  inline const heap_type & get_derivations() const { return derivations; }
  
  inline const packed_edge_probability_with_index& get(unsigned idx) const {return derivations[idx];}
  inline packed_edge_probability& get(unsigned idx) { return derivations[idx]; }
  
  template<class TDaughter>
  inline void update_best(const TDaughter& dtr);
  inline void finalize_best();
  
  inline void find_succ(packed_edge_probability_with_index& pep, bool licence_unaries);
  void extend_derivation(unsigned, bool) ;
  
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
};

inline std::ostream& operator<<(std::ostream& out, const MinDivBest & best)
{
  return out << "(Best: " << &best << ")";
}

class MinDivProbabilityKB
{

public:

  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::PEdge PEdge;
  typedef typename MinDivKBTypes::UEdge UEdge;
  typedef typename MinDivKBTypes::LBEdge LBEdge;
  typedef typename MinDivKBTypes::Cell Cell;
  typedef typename MinDivKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MinDivKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MinDivKBTypes::LexicalDaughter LexicalDaughter;

private:
  friend class ParserCKYAllMinDivKB;

  double inside_q;
  double outside_q;
  double inside_p;
  double outside_p;

  static double log_normalisation_factor;
  static double normalisation_factor;
  static double log_normalisation_factor_q;
  static double     normalisation_factor_q;


  MinDivProbabilityKB();
public:
  ~MinDivProbabilityKB() {}

  inline static void   set_normalisation_factor(double lnf) {normalisation_factor = lnf;};
  inline static double get_normalisation_factor() {return normalisation_factor;};
  inline static void   set_log_normalisation_factor(double lnf) {log_normalisation_factor = lnf;};
  inline static double get_log_normalisation_factor() {return log_normalisation_factor;};
  inline static void   set_log_normalisation_factor_q(double lnf) {log_normalisation_factor_q = lnf;};
  inline static double get_log_normalisation_factor_q() {return log_normalisation_factor_q;};
  inline static void   set_normalisation_factor_q(double lnf) {normalisation_factor_q = lnf;};
  inline static double get_normalisation_factor_q() {return normalisation_factor_q;};
//   inline static void   set_normalisation_factor(double lnf) {normalisation_factor = lnf;};
//   inline static double get_normalisation_factor() {return normalisation_factor;};
  inline double & get_inside_q() { return inside_q; }
  inline const double & get_inside_q() const { return inside_q; }
  inline double & get_outside_q() { return outside_q; }
  inline const double & get_outside_q() const { return outside_q; }
  void set_outside_q(double prob) {outside_q = prob;}
  

  inline void reinit_inside_outside(double val);
  inline void reinit_q_inside_outside(double val);

  inline void update_inside_q_lexical(const LexicalDaughter& dtr);
  inline void update_inside_q_unary(const UnaryDaughter& dtr);
  inline void update_inside_q_binary(const BinaryDaughter& dtr);
  
  inline void update_outside_q_lexical(LexicalDaughter& dtr);
  inline void update_outside_q_unary(UnaryDaughter& dtr);
  inline void update_outside_q_binary(BinaryDaughter& dtr);
  
  inline void update_q_lexical(LexicalDaughter& dtr);
  inline void update_q_unary(UnaryDaughter& dtr);
  inline void update_q_binary(BinaryDaughter& dtr);

  
  
  public:
    inline std::ostream& operator>>(std::ostream& out) const;
};

inline std::ostream& operator<<(std::ostream& out, const MinDivProbabilityKB & prob)
{
  return out << "(MinDivProb: " << &prob << ")";
}


#include "MinDivDaughters.h"


class ParserCKYAllMinDivKB : public ParserCKYAll_Impl<MinDivKBTypes>
{
  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::Cell Cell;
  typedef typename MinDivKBTypes::Best Best;
  typedef typename MinDivKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MinDivKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MinDivKBTypes::LexicalDaughter LexicalDaughter;

private:
  unsigned k;

public:
  ParserCKYAllMinDivKB(std::vector<AGrammar*>& cgs,
                       const std::vector<double>& p, double b_t,
                       const annot_descendants_type& annot_descendants_,
                       bool accurate_, unsigned min_beam, int stubborn, unsigned k_);

  inline ~ParserCKYAllMinDivKB() {};

  void extract_solution();

private:
  /** also computes marginals for each daughter */
  inline virtual void compute_outside_probabilities();
  
  /* computes inside-outside on q */
  inline void compute_inside_q_probabilities();
  inline void compute_outside_q_probabilities();
  inline void compute_inside_outside_q_probabilities();
  inline double get_sentence_probability_q() const;

  /* computes q as marginal(p) / (inside(q)*outside(q)) */
  inline void update_q();
  
  /* filling edge probability structures with "best" pointers */
  inline void fill_bests();
  inline void initialise_candidates();

  inline void extend_all_derivations();
};


template<>
std::ostream & operator<<(std::ostream &, const BasePackedEdge<MinDivKBTypes> &);




#endif /* _PARSERCKYALLMINDIVKB_H_ */
