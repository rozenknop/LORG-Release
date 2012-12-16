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

template<class Types> class UPackedEdge;
template<class Types> class LBPackedEdge;
template<class Types> class BasePackedEdge;

template<class Types>
class BasePackedEdge
{
public:
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Best Best;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef BasePackedEdge<Types> PEdge;
  typedef UPackedEdge<Types> UEdge;
  typedef LBPackedEdge<Types> LBEdge;
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

protected :
  /**
   * attributes
   */
  Best best;
  AnnotationInfo annotations;  ///< probabilities
  bool open;
  ProbaModel proba;


private:

  /**
     \brief Constructors are forbidden
  */
  BasePackedEdge() {}
  BasePackedEdge(const BasePackedEdge<Types> &) {}




public:

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
  friend std::ostream& operator<<(std::ostream& out, const BasePackedEdge<OPEP>& edge);

  void dump(std::ostream & out) const;

  /**
     \brief return an AnnotationInfo (inside/outside probs +scaling values)
   */
  AnnotationInfo& get_annotations();
  const AnnotationInfo& get_annotations() const;



  /**
     \brief get the structure holding the "best calculation for the prob. model" whatever it means
   */
  inline const ProbaModel& get_prob_model() const;
  inline ProbaModel& get_prob_model();

  inline bool valid_prob_at(unsigned i) const;

  inline bool is_closed() const { return not open; }
  inline bool is_opened() const { return open; }

  inline Best& get_best();
  inline const Best& get_best() const;
  inline void extend_derivation(unsigned i, bool licence_unaries);
  inline bool has_solution(unsigned i) const ;

  PtbPsTree * to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const;

protected:
  void to_ptbpstree(PtbPsTree& tree, PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
                    bool append_annot, bool outpu_forms) const;

public:
  
  void process(function<void(ProbaModel &, const AnnotationInfo &)> f) {f(get_prob_model(), get_annotations());}
  void process(function<void(ProbaModel &)> f) {f(get_prob_model());}
  void process(function<void(Best &, const AnnotationInfo &)> f) {f(get_best(), get_annotations());}
  void process(function<void(AnnotationInfo &)> f) {f(get_annotations());}
  void process(function<void(Best &)> f) {f(get_best());}
  void process(function<void(BasePackedEdge &)> f) { f(*this); }
};


template<class Types>
class UPackedEdge : public BasePackedEdge<Types>
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

protected:
  uvector unary_daughters;     ///< set of possible daughters
  static PathMatrix unary_chains;

public:
  /**
   * \brief replace rules with grammar #i
   */
  void replace_rule_probabilities(unsigned i);

  inline const uvector& get_unary_daughters() const;
  inline uvector& get_unary_daughters();

 /**
     \brief get one particular daughter of the edge
  */
  inline UnaryDaughter& get_unary_daughter(unsigned i);
  inline const UnaryDaughter& get_unary_daughter(unsigned i) const;

  /**
   *   \brief build and add a daughter
   */
  inline void add_daughters(LBPackedEdge<Types> & left, const UnaryRule* rule);

  /**
   * \brief set unary chains
   */
  static void set_unary_chains(const PathMatrix& pathmatrix);
  static const PathMatrix& get_unary_chains();

  void dump(std::ostream & out) const;


  inline bool no_daughters() { return unary_daughters.empty(); } 
  void close() { this->open=false; UPackedEdge::~UPackedEdge(); }
};



template<class Types>
class LBPackedEdge : public BasePackedEdge<Types>
{
public:
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::PEdge PEdge;
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

protected:
  bvector binary_daughters;    ///< set of possible daughters
  lvector lexical_daughters;   ///< a possible lexical daughter

public:
  void reserve_binary_daughters(int size) 
  {
    binary_daughters.reserve(size);
  }

  /**
   * \brief replace rules with grammar #i
   */
  void replace_rule_probabilities(unsigned i);

  /**
   *   \brief get the daughters of the edge
   */
  const bvector& get_binary_daughters() const;
  bvector& get_binary_daughters();

  const lvector& get_lexical_daughters() const;
  lvector& get_lexical_daughters();

 /**
     \brief get one particular daughter of the edge
  */
  BinaryDaughter& get_binary_daughter(unsigned i);
  const BinaryDaughter& get_binary_daughter(unsigned i) const;

  LexicalDaughter& get_lexical_daughter(unsigned i);
  const LexicalDaughter& get_lexical_daughter(unsigned i) const;

  /**
     \brief is there a lexical daughter ?
  */
  bool get_lex() const;

  /**
   *   \brief build and add a daughter (binary and lexical versions)
   */
  inline void add_daughters(PEdge & left, PEdge & right, const BinaryRule* rule);
  inline void add_daughters(const LexicalRule* rule, const Word* w);

  void clean_invalidated_binaries();

  bool no_daughters() { return binary_daughters.empty() and lexical_daughters.empty(); } 
  void close() { this->open=false; LBPackedEdge<Types>::~LBPackedEdge(); }


  PtbPsTree * to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const
  { 
    return BasePackedEdge<Types>::to_ptbpstree(lhs,ith_deriv,append_annot,output_forms);
  }

private:
  void to_ptbpstree(PtbPsTree& tree, PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
                    bool append_annot, bool outpu_forms) const;
                    
  friend class BasePackedEdge<Types>;
  
public:
    void dump(std::ostream & out) const;
};



template<class Types>
class PackedEdge
{
public:
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef BasePackedEdge<Types> PEdge;
  typedef UPackedEdge<Types> UEdge;
  typedef LBPackedEdge<Types> LBEdge;
  typedef typename Types::Best Best;

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
  UEdge u;
  LBEdge lb;
  
public:

  LBEdge & lbedge() { return lb; }
  UEdge  &  uedge() { return u; }
  const LBEdge & lbedge() const { return lb; }
  const UEdge  &  uedge() const { return u; }
  
  
  inline bool has_solution(unsigned i) const { return lb.has_solution() or u.has_solution(); }
  inline bool is_closed() const { return lb.is_closed() and u.is_closed(); }

  
  void process(function<void(PEdge &)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  void process(function<void(Best&)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  void process(function<void(ProbaModel&)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  void process(function<void(Best&, const AnnotationInfo&)> f) { 
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  inline void process(function<void(AnnotationInfo&)> f) { 
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  
  void process(function<void(UEdge &)> f) { if (u.is_opened()) f(u); }
  void process(function<void(const UnaryDaughter &)> f) const { if (u.is_opened()) for(const auto& d: u.get_unary_daughters()) f(d); }
  void process(function<void(UEdge &, UnaryDaughter &)> f) { if (u.is_opened()) for(auto& d: u.get_unary_daughters()) f(u, d); }
  void process(function<void(const UnaryDaughter &, AnnotationInfo &)> f) { if (u.is_opened()) for(const auto& d: u.get_unary_daughters()) f(d, u.get_annotations()); }
  void process(function<void(UnaryDaughter &, AnnotationInfo &)> f) { if (u.is_opened()) for(auto& d: u.get_unary_daughters()) f(d, u.get_annotations()); }
  void process(function<void(Best &, UEdge &, const UnaryDaughter &)> f) {if (u.is_opened()) for(const auto& d: u.get_unary_daughters()) f(u.get_best(), u, d);}
  void process(function<void(ProbaModel &, const UnaryDaughter &)> f) {if (u.is_opened()) for(const auto& d: u.get_unary_daughters()) f(u.get_prob_model(), d);}
  void process(function<void(ProbaModel &, UnaryDaughter &)> f) {if (u.is_opened()) for(auto& d: u.get_unary_daughters()) f(u.get_prob_model(), d);}
  void process(function<void(Best &, const UnaryDaughter &)> f) {if (u.is_opened()) for(const auto& d: u.get_unary_daughters()) f(u.get_best(), d);}
  void process(function<void(Best &, UnaryDaughter &)> f) {if (u.is_opened()) for(auto& d: u.get_unary_daughters()) f(u.get_best(), d);}

  void process(function<void(LBEdge &)> f) { if (lb.is_opened()) f(lb); }

  void process(function<void(const LexicalDaughter &)> f) const {if (lb.is_opened()) for(const auto& d: lb.get_lexical_daughters()) f(d);}
  void process(function<void(const BinaryDaughter &)> f) const { if (lb.is_opened()) for(const auto& d: lb.get_binary_daughters()) f(d); }

  void process(function<void(LBEdge &, LexicalDaughter &)> f) {if (lb.is_opened()) for(auto& d: lb.get_lexical_daughters()) f(lb, d);}
  void process(function<void(LBEdge &, BinaryDaughter &)> f) { if (lb.is_opened()) for(auto& d: lb.get_binary_daughters()) f(lb, d); }

  void process(function<void(const LexicalDaughter &, AnnotationInfo &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_lexical_daughters()) f(d, lb.get_annotations());}
  void process(function<void(const BinaryDaughter &, AnnotationInfo &)> f) { if (lb.is_opened()) for(const auto& d: lb.get_binary_daughters()) f(d, lb.get_annotations()); }

  void process(function<void(LexicalDaughter &, AnnotationInfo &)> f) {if (lb.is_opened()) for(auto& d: lb.get_lexical_daughters()) f(d, lb.get_annotations());}
  void process(function<void(BinaryDaughter &, AnnotationInfo &)> f) { if (lb.is_opened()) for(auto& d: lb.get_binary_daughters()) f(d, lb.get_annotations()); }

  void process(function<void(Best &, LBEdge &, const LexicalDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_lexical_daughters()) f(lb.get_best(), lb, d);}
  void process(function<void(Best &, LBEdge &, const BinaryDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_binary_daughters()) f(lb.get_best(), lb, d);}

  void process(function<void(ProbaModel &, const LexicalDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_lexical_daughters()) f(lb.get_prob_model(), d);}
  void process(function<void(ProbaModel &, const BinaryDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_binary_daughters()) f(lb.get_prob_model(), d);}
  void process(function<void(ProbaModel &, LexicalDaughter &)> f) {if (lb.is_opened()) for(auto& d: lb.get_lexical_daughters()) f(lb.get_prob_model(), d);}
  void process(function<void(ProbaModel &, BinaryDaughter &)> f) {if (lb.is_opened()) for(auto& d: lb.get_binary_daughters()) f(lb.get_prob_model(), d);}

  void process(function<void(Best &, const LexicalDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_lexical_daughters()) f(lb.get_best(), d);}
  void process(function<void(Best &, const BinaryDaughter &)> f) {if (lb.is_opened()) for(const auto& d: lb.get_binary_daughters()) f(lb.get_best(), d);}
  void process(function<void(Best &, LexicalDaughter &)> f) {if (lb.is_opened()) for(auto& d: lb.get_lexical_daughters()) f(lb.get_best(), d);}
  void process(function<void(Best &, BinaryDaughter &)> f) {if (lb.is_opened()) for(auto& d: lb.get_binary_daughters()) f(lb.get_best(), d);}

  void apply() const {}
  template<typename Function, typename... OtherFunctions>
  void apply(Function&& f, OtherFunctions&&... o) {process(toFunc(f));apply(o...);}
  template<typename Function, typename... OtherFunctions>
  void apply(Function&& f, OtherFunctions&&... o) const {process(toFunc(f));apply(o...);}

  template<class OPEP>
  friend std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge);
  void dump(std::ostream & out) const;
};

template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const BasePackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const UPackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const LBPackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge) {edge.dump(out); return out; }


#endif /*PACKEDEDGE_H_*/
