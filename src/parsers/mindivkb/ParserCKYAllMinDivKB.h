// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMINDIVKB_H_
#define _PARSERCKYALLMINDIVKB_H_

#include "parsers/ParserCKYAll.h"
#include "emptystruct.h"
#include "MinDivTypes.h"



class MinDivProbabilityKB
{

public:

  typedef std::vector<packed_edge_probability_with_index> heap_type;

  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::Cell Cell;
  typedef typename MinDivKBTypes::UnaryDaughter UnaryDaughter;
  typedef typename MinDivKBTypes::BinaryDaughter BinaryDaughter;
  typedef typename MinDivKBTypes::LexicalDaughter LexicalDaughter;

private:
  friend class ParserCKYAllMinDivKB;

  heap_type candidates;
  heap_type derivations;

  double inside_prob;
  double outside_prob;
  double inside_unary_temp ;
  double outside_unary_temp ;

  static double log_normalisation_factor;
  static double normalisation_factor;
  static unsigned size;

public:

  MinDivProbabilityKB() :  candidates(), derivations() {candidates.reserve(50);};
  ~MinDivProbabilityKB() {};

  inline static void set_size(unsigned k) {size = k;}

  inline static void   set_log_normalisation_factor(double lnf) {log_normalisation_factor = lnf;};
  inline static double get_log_normalisation_factor() {return log_normalisation_factor;};
  inline static void   set_normalisation_factor(double lnf) {normalisation_factor = lnf;};
  inline static double get_normalisation_factor() {return normalisation_factor;};
  inline double & get_inside_prob() { return inside_prob; }
  inline const double & get_inside_prob() const { return inside_prob; }
  inline double & get_outside_prob() { return outside_prob; }
  void set_outside_prob(double prob) {outside_prob = prob;}
  inline const double & get_outside_prob() const { return outside_prob; }
  inline double & get_inside_unary_temp() { return inside_unary_temp; }
  inline const double & get_inside_unary_temp() const { return inside_unary_temp; }
  inline double & get_outside_unary_temp() { return outside_unary_temp; }
  inline const double & get_outside_unary_temp() const { return outside_unary_temp; }
  
  inline const packed_edge_probability_with_index& get(unsigned idx) const {return derivations[idx];}
  inline packed_edge_probability& get(unsigned idx) { return derivations[idx]; }

  template<class TDaughter>
  inline void update_best(const TDaughter& dtr);
  inline void finalize_best();
  
  inline void find_succ(packed_edge_probability_with_index& pep, bool licence_unaries);
  void extend_derivation(unsigned, bool) ;

  inline unsigned n_deriv() const {return derivations.size();};

  inline bool has_solution(unsigned i) const {return i <derivations.size();}

  inline void reinit_inside_outside(double val);

  inline void update_inside_lexical(const LexicalDaughter& dtr);
  inline void prepare_inside_unary();
  inline void update_inside_unary(const UnaryDaughter& dtr);
  inline void adjust_inside_unary();
  inline void update_inside_binary(const BinaryDaughter& dtr);
  
  inline void update_outside_lexical(const LexicalDaughter& dtr);
  inline void prepare_outside_unary();
  inline void update_outside_unary(const UnaryDaughter& dtr);
  inline void adjust_outside_unary();
  inline void update_outside_binary(const BinaryDaughter& dtr);
  
  inline void update_q_lexical(LexicalDaughter& dtr);
  inline void update_q_unary(UnaryDaughter& dtr);
  inline void update_q_binary(BinaryDaughter& dtr);

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



#include "MinDivDaughters.h"


class ParserCKYAllMinDivKB : public ParserCKYAll_Impl<MinDivKBTypes>
{
  typedef typename MinDivKBTypes::Edge Edge;
  typedef typename MinDivKBTypes::Cell Cell;
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

  ~ParserCKYAllMinDivKB() {};

  void extract_solution();

private:
  /** also computes marginals for each daughter */
  virtual void compute_outside_probabilities();
  
  /* computes inside-outside on q */
  void compute_inside_q_probabilities();
  void compute_outside_q_probabilities();
  void compute_inside_outside_q_probabilities();

  /* computes q as marginal(p) / (inside(q)*outside(q)) */
  void update_q();
  
  /* filling edge probability structures with "best" pointers */
  void fill_bests();
  void initialise_candidates();

  void extend_all_derivations();
};

#endif /* _PARSERCKYALLMINDIVKB_H_ */
