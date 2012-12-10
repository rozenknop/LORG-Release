// -*- mode: c++ -*-
#ifndef PACKEDEDGE_H_
#define PACKEDEDGE_H_

#include "PackedEdgeDaughters.h"

#include "grammars/AnnotatedLabelsInfo.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include <unordered_map>

#include "AnnotationInfo.h"

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"

#include "utils/PtbPsTree.h"

#include "PCKYAllCell.h"


#include "utils/SymbolTable.h"
#include "PackedEdgeProbability.h"

#include "utils/lorg_functional.h"
#include "utils/hash_impl.h"
#include <utils/tick_count.h>


typedef std::pair< int, unsigned> asymb;
typedef std::unordered_map<asymb, std::unordered_map< asymb, asymb> > PathMatrix;


/**
  \class PackedEdge
  \brief represents an edge in a chart
*/



template<class Types>
class PackedEdge
{
public:
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
//   typedef PackedEdgeDaughters Daughters;
  typedef typename Types::BRule BinaryRule;
  typedef typename Types::URule UnaryRule;
  typedef typename Types::LRule LexicalRule;
  
  typedef typename Types::BinaryDaughter  BinaryDaughter;
  typedef typename Types::UnaryDaughter   UnaryDaughter;
  typedef typename Types::LexicalDaughter LexicalDaughter;

  typedef std::vector<BinaryDaughter> bvector;
  typedef std::vector<UnaryDaughter>  uvector;
  typedef std::vector<LexicalDaughter>  lvector;

  typedef typename bvector::iterator biterator;
  typedef typename uvector::iterator  uiterator;
  typedef typename lvector::iterator  literator;

  typedef typename bvector::const_iterator const_biterator;
  typedef typename uvector::const_iterator  const_uiterator;
  typedef typename lvector::const_iterator  const_literator;

private:

  /**
     \brief Forbidden constructors
  */
  PackedEdge() {}
  PackedEdge(const PackedEdge<Types> &o) {}

public:
  void reserve_binary_daughters(int size) 
  {
    binary_daughters.reserve(size);
  }

  /**
     \brief Destructor
  */
  //   ~PackedEdge()  {assert(false);}


  /**
     \brief get the daughters of the edge
  */
  const bvector& get_binary_daughters() const;
  bvector& get_binary_daughters();

  const uvector& get_unary_daughters() const;
  uvector& get_unary_daughters();

  const lvector& get_lexical_daughters() const;
  lvector& get_lexical_daughters();



 /**
     \brief get one particular daughter of the edge
  */
  BinaryDaughter& get_binary_daughter(unsigned i);
  const BinaryDaughter& get_binary_daughter(unsigned i) const;

  UnaryDaughter& get_unary_daughter(unsigned i);
  const UnaryDaughter& get_unary_daughter(unsigned i) const;

  LexicalDaughter& get_lexical_daughter(unsigned i);
  const LexicalDaughter& get_lexical_daughter(unsigned i) const;

  /**
     \brief is the daughter lexical ?
  */
  bool get_lex() const;

  /**
     \brief resize the vector of annotations
     \param size the new size
   */
  void local_resize_annotations(unsigned size);

  /**
     \brief Output operator
     \param out the ostream to write on
     \param edge the edge object to write
     \return the used ostream
  */
  template<class OPEP>
  friend std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge);



  /**
     \brief return an AnnotationInfo (inside/outside probs +scaling values)
   */
  AnnotationInfo& get_annotations();
  const AnnotationInfo& get_annotations() const;



  /**
   *  \brief before and after computing the inside probability of a packed edge
   */
  void prepare_inside_probability();
  void adjust_inside_probability();


  /**
   * \brief before and after computing the outside probability of a packed edge
   */
  void prepare_outside_probability();
  void adjust_outside_probability();

  /**
     \brief build and add a daughter (binary, unary and lexical versions)
   */
  void add_daughters(Edge & left,
                     Edge & right, const BinaryRule* rule);
  void add_daughters(Edge & left, const UnaryRule* rule);
  void add_daughters(const LexicalRule* rule, const Word* w);

  /**
     \brief get the structure holding the "best calculation for the prob. model" whatever it means
   */
  const ProbaModel& get_prob_model() const;
  ProbaModel& get_prob_model();


  /**
     \brief set unary chains
   */
  static void set_unary_chains(const PathMatrix& pathmatrix);
  static const PathMatrix& get_unary_chains();


  void replace_rule_probabilities(unsigned i);


  void extend_derivation(unsigned i, bool licence_unaries);

  bool valid_prob_at(unsigned i) const;

  void clean_invalidated_binaries();

  PtbPsTree * to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const;

  bool has_solution(unsigned i) const ;

  bool no_daughters() { return binary_daughters.empty() and unary_daughters.empty() and lexical_daughters.empty(); } 
  bool is_closed() const { return not open; }
  void close() { open=false; Edge::~Edge(); }
  
protected :
  bvector binary_daughters;    ///< set of possible daughters
  uvector unary_daughters;     ///< set of possible daughters
  lvector lexical_daughters;   ///< a possible lexical daughter
  AnnotationInfo annotations;  ///< probabilities

  static PathMatrix unary_chains;

  bool open;

private:
  ProbaModel best;

  void to_ptbpstree(PtbPsTree& tree, PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
                    bool append_annot, bool outpu_forms) const;

public:
  void process(function<void(const LexicalDaughter &)> f) const {for(const auto& d: get_lexical_daughters()) f(d);}
  void process(function<void(const UnaryDaughter &)> f) const { for(const auto& d: get_unary_daughters()) f(d); }
  void process(function<void(const BinaryDaughter &)> f) const { for(const auto& d: get_binary_daughters()) f(d); }

  void process(function<void(Edge &, LexicalDaughter &)> f) {for(auto& d: get_lexical_daughters()) f(*this, d);}
  void process(function<void(Edge &, UnaryDaughter &)> f) { for(auto& d: get_unary_daughters()) f(*this, d); }
  void process(function<void(Edge &, BinaryDaughter &)> f) { for(auto& d: get_binary_daughters()) f(*this, d); }

  void process(function<void(const LexicalDaughter &, AnnotationInfo &)> f) {for(const auto& d: get_lexical_daughters()) f(d, get_annotations());}
  void process(function<void(const UnaryDaughter &, AnnotationInfo &)> f) { for(const auto& d: get_unary_daughters()) f(d, get_annotations()); }
  void process(function<void(const BinaryDaughter &, AnnotationInfo &)> f) { for(const auto& d: get_binary_daughters()) f(d, get_annotations()); }

  void process(function<void(LexicalDaughter &, AnnotationInfo &)> f) {for(auto& d: get_lexical_daughters()) f(d, get_annotations());}
  void process(function<void(UnaryDaughter &, AnnotationInfo &)> f) { for(auto& d: get_unary_daughters()) f(d, get_annotations()); }
  void process(function<void(BinaryDaughter &, AnnotationInfo &)> f) { for(auto& d: get_binary_daughters()) f(d, get_annotations()); }
  
  void process(function<void(ProbaModel &, const AnnotationInfo &)> f) {f(get_prob_model(), get_annotations());}

  void process(function<void(ProbaModel &, Edge &, const LexicalDaughter &)> f) {for(const auto& d: get_lexical_daughters()) f(get_prob_model(), *this, d);}
  void process(function<void(ProbaModel &, Edge &, const UnaryDaughter &)> f) {for(const auto& d: get_unary_daughters()) f(get_prob_model(), *this, d);}
  void process(function<void(ProbaModel &, Edge &, const BinaryDaughter &)> f) {for(const auto& d: get_binary_daughters()) f(get_prob_model(), *this, d);}

  void process(function<void(ProbaModel &, const LexicalDaughter &)> f) {for(const auto& d: get_lexical_daughters()) f(get_prob_model(), d);}
  void process(function<void(ProbaModel &, const UnaryDaughter &)> f) {for(const auto& d: get_unary_daughters()) f(get_prob_model(), d);}
  void process(function<void(ProbaModel &, const BinaryDaughter &)> f) {for(const auto& d: get_binary_daughters()) f(get_prob_model(), d);}

  void process(function<void(ProbaModel &, LexicalDaughter &)> f) {for(auto& d: get_lexical_daughters()) f(get_prob_model(), d);}
  void process(function<void(ProbaModel &, UnaryDaughter &)> f) {for(auto& d: get_unary_daughters()) f(get_prob_model(), d);}
  void process(function<void(ProbaModel &, BinaryDaughter &)> f) {for(auto& d: get_binary_daughters()) f(get_prob_model(), d);}
  
  void process(function<void(ProbaModel &)> f) {f(get_prob_model());}

  void process(function<void(Edge &)> f) { f(*this); }
  
  template<typename Function, typename... OtherFunctions>
  void apply(Function&& f, OtherFunctions&&... o) {process(toFunc(f));apply(o...);}
  void apply() const {}
  template<typename Function, typename... OtherFunctions>
  void apply(Function&& f, OtherFunctions&&... o) const {process(toFunc(f));apply(o...);}
};

#endif /*PACKEDEDGE_H_*/
