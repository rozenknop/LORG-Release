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
  \class AnnotatedEdge
  \brief represents inside and outside annotations of an edge in a chart
*/

template<class Types> class UPackedEdge;
template<class Types> class LBPackedEdge;
template<class Types> class BasePackedEdge;
template<class Types> class PackedEdge;

template<class Types>
class AnnotatedEdge
{
public:
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Best Best;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef AnnotatedEdge<Types> AEdge;
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
  AnnotationInfo annotations;  ///< probabilities
  bool open;
  friend class PackedEdge<Types>;

private:

  /**
     \brief Constructors are forbidden
  */
  AnnotatedEdge() {}
  AnnotatedEdge(const AnnotatedEdge<Types> &) {}




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
  friend std::ostream& operator<<(std::ostream& out, const AnnotatedEdge<OPEP>& edge);

  void dump(std::ostream & out) const;

  /**
     \brief return an AnnotationInfo (inside/outside probs +scaling values)
   */
  AnnotationInfo& get_annotations();
  const AnnotationInfo& get_annotations() const;
  inline bool valid_prob_at(unsigned i) const;

  inline bool is_closed() const { return not open; }
  inline bool is_opened() const { return open; }


public:
  
  inline void process(function<void(AnnotationInfo &)> f) {f(get_annotations());}
};


/**
  \class BasePackedEdge
  \brief represents an edge in a chart
*/

template<class Types> class UPackedEdge;
template<class Types> class LBPackedEdge;
template<class Types> class BasePackedEdge;

template<class Types>
class BasePackedEdge : public AnnotatedEdge<Types>
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
  ProbaModel proba;


private:

  /**
     \brief Constructors are forbidden
  */
  BasePackedEdge() {}
  BasePackedEdge(const BasePackedEdge<Types> &) {}




public:

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
     \brief get the structure holding the "best calculation for the prob. model" whatever it means
   */
  inline const ProbaModel& get_prob_model() const;
  inline ProbaModel& get_prob_model();

  inline Best& get_best();
  inline const Best& get_best() const;
  inline void extend_derivation(unsigned i, bool licence_unaries);
  inline bool has_solution(unsigned i) const ;

  PtbPsTree * to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const;

protected:
  void to_ptbpstree(PtbPsTree& tree, PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
                    bool append_annot, bool outpu_forms) const;

public:
  
  inline void process(function<void(ProbaModel &, const AnnotationInfo &)> f) {f(get_prob_model(), AnnotatedEdge<Types>::get_annotations());}
  inline void process(function<void(ProbaModel &)> f) {f(get_prob_model());}
  inline void process(function<void(Best &, const AnnotationInfo &)> f) {f(get_best(), AnnotatedEdge<Types>::get_annotations());}
  inline void process(function<void(Best &)> f) {f(get_best());}
  inline void process(function<void(BasePackedEdge &)> f) { f(*this); }
  inline void process(function<void(AnnotationInfo &)> f) {AnnotatedEdge<Types>::process(f);}
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
  friend class PackedEdge<Types>;
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
   * \brief set unary chains
   */
  static void set_unary_chains(const PathMatrix& pathmatrix);
  static const PathMatrix& get_unary_chains();

  void dump(std::ostream & out) const;


  inline bool no_daughters() { return unary_daughters.empty(); } 
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
  friend class PackedEdge<Types>;

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

  void clean_invalidated_binaries();

  bool no_daughters() { return binary_daughters.empty() and lexical_daughters.empty(); } 


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
class PackedEdge : public AnnotatedEdge<Types>
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
  
  /**
   *   \brief build and add a daughter (lexical version)
   */
  inline void add_daughters(Edge & left, const UnaryRule* rule);
  /**
   *   \brief build and add a daughter (lexical and binary versions)
   */
  inline void add_daughters(const LexicalRule* rule, const Word* w);
  inline void add_daughters(Edge & left, Edge & right, const BinaryRule* rule);

  /**
   * \brief close the edge, releasing memory
   */
  inline void close() { 
    this->open = lb.open = u.open = false;
    PackedEdge::~PackedEdge();
  }
  /**
   * \brief close the unary edge, and the PackedEdge if lexical/binary part already closed
   */
  inline void close_u() { 
    if (lb.is_closed()) 
      close();
    else 
      u.open = false;
  }
  /**
   * \brief close the lexical/binary edge, and the PackedEdge if unary part already closed
   */
  inline void close_lb() { 
    if (u.is_closed()) 
      close();
    else 
      lb.open = false;
  }
  
  inline void reset_probabilities() {
    this->get_annotations().reset_probabilities();
    if (lb.is_opened()) lb.get_annotations().reset_probabilities();
    if (u.is_opened()) u.get_annotations().reset_probabilities();
  }
  inline void resize_annotations(unsigned size) {
    this->get_annotations().reset_probabilities();
    this->get_annotations().resize(size);
    if (lb.is_opened()) {
      lb.get_annotations().reset_probabilities();
      lb.get_annotations().resize(size);
    }
    if (u.is_opened()) {
      u.get_annotations().reset_probabilities();
      u.get_annotations().resize(size);
    }
  }
  /**
   * inside computation
   */
  inline void update_merged_inside_annotations_from_lex() {
    for(const auto& d: lb.get_lexical_daughters())
      d.update_inside_annotations(this->annotations);
  }
  inline void update_merged_inside_annotations_from_bin() {
    for(const auto& d: lb.get_binary_daughters()) {
      d.update_inside_annotations(this->annotations);
    }
  }
  inline void update_unary_inside_annotations() {
    for(const auto& d: u.get_unary_daughters()) {
      d.update_inside_annotations(u.annotations);
    }
  }
  inline void update_binary_inside_annotations() {
    for(const auto& d: lb.get_binary_daughters()) {
      d.update_inside_annotations(lb.annotations);
    }
  }
  inline void update_lexical_inside_annotations() {
    for(const auto& d: lb.get_lexical_daughters()) {
      d.update_inside_annotations(lb.annotations);
    }
  }
  inline void add_unary_insides_to_merged() {
    for(unsigned i=0; i<this->annotations.get_size(); ++i)
      if (u.annotations.inside_probabilities.array[i]!=LorgConstants::NullProba)
        this->annotations.inside_probabilities.array[i] += u.annotations.inside_probabilities.array[i];
  }

  /**
   * outside computation
   */
  inline void copy_merged_outsides_to_unary() {
    for(unsigned i=0; i<this->annotations.get_size(); ++i)
      u.annotations.outside_probabilities.array[i] = this->annotations.outside_probabilities.array[i];
  }
  inline void update_unary_outside_annotations() {
    for(const auto& d: u.get_unary_daughters()) {
      d.update_outside_annotations(u.annotations);
    }
  }
  inline void update_binary_outside_annotations() {
    for(const auto& d: lb.get_binary_daughters()) {
      d.update_outside_annotations(lb.annotations);
    }
  }
  inline void update_lexical_outside_annotations() {
    for(const auto& d: lb.get_lexical_daughters()) {
      d.update_outside_annotations(lb.annotations);
    }
  }
  inline void update_binary_outside_annotations_from_merged() {
    for(const auto& d: lb.get_binary_daughters()) {
      d.update_outside_annotations(this->annotations);
    }
  }
  


//   inline void add_from_lb_insides() {
//     for(unsigned i=0; i<this->annotations.get_size(); ++i)
//       this->annotations.inside_probabilities.array[i] += lb.annotations.inside_probabilities.array[i];
//   }
//   inline void add_from_lb_outsides() {
//     for(unsigned i=0; i<this->annotations.get_size(); ++i)
//       this->annotations.outside_probabilities.array[i] += lb.annotations.outside_probabilities.array[i];
//   }
//   inline void add_to_lb_outsides() {
//     for(unsigned i=0; i<this->annotations.get_size(); ++i)
//       lb.annotations.outside_probabilities.array[i] += this->annotations.outside_probabilities.array[i];
//   }

  
  inline void process(function<void(PEdge &)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  inline void process(function<void(Best&)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  inline void process(function<void(ProbaModel&)> f) {
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  inline void process(function<void(Best&, const AnnotationInfo&)> f) { 
    if (lb.is_opened()) lb.process(f);
    if (u.is_opened()) u.process(f); }
  
  
  inline void process(function<void(Edge &)> f) { f(*this); }

  inline void process(function<void(UEdge &)> f) { f(u); }
  inline void process(function<void(const UnaryDaughter &)> f)                   const { for(const auto& d: u.get_unary_daughters()) f(d); }
  inline void process(function<void(UEdge &, UnaryDaughter &)> f)                      { for(auto& d: u.get_unary_daughters()) f(u, d); }
  inline void process(function<void(const UnaryDaughter &, AnnotationInfo &)> f)       { for(const auto& d: u.get_unary_daughters()) f(d, u.get_annotations()); }
  inline void process(function<void(UnaryDaughter &, AnnotationInfo &)> f)             { for(auto& d: u.get_unary_daughters()) f(d, u.get_annotations()); }
  inline void process(function<void(Best &, UEdge &, const UnaryDaughter &)> f)        { for(const auto& d: u.get_unary_daughters()) f(u.get_best(), u, d);}
  inline void process(function<void(ProbaModel &, const UnaryDaughter &)> f)           { for(const auto& d: u.get_unary_daughters()) f(u.get_prob_model(), d);}
  inline void process(function<void(ProbaModel &, UnaryDaughter &)> f)                 { for(auto& d: u.get_unary_daughters()) f(u.get_prob_model(), d);}
  inline void process(function<void(Best &, const UnaryDaughter &)> f)                 { for(const auto& d: u.get_unary_daughters()) f(u.get_best(), d);}
  inline void process(function<void(Best &, UnaryDaughter &)> f)                       { for(auto& d: u.get_unary_daughters()) f(u.get_best(), d);}

  inline void process(function<void(LBEdge &)> f) { f(lb); }

  inline void process(function<void(const LexicalDaughter &)> f) const {for(const auto& d: lb.get_lexical_daughters()) f(d);}
  inline void process(function<void(const BinaryDaughter &)> f) const { for(const auto& d: lb.get_binary_daughters()) f(d); }

  inline void process(function<void(LBEdge &, LexicalDaughter &)> f) {for(auto& d: lb.get_lexical_daughters()) f(lb, d);}
  inline void process(function<void(LBEdge &, BinaryDaughter &)> f) { for(auto& d: lb.get_binary_daughters()) f(lb, d); }

  inline void process(function<void(const LexicalDaughter &, AnnotationInfo &)> f) {for(const auto& d: lb.get_lexical_daughters()) f(d, lb.get_annotations());}
  inline void process(function<void(const BinaryDaughter &, AnnotationInfo &)> f) { for(const auto& d: lb.get_binary_daughters()) f(d, lb.get_annotations()); }

  inline void process(function<void(LexicalDaughter &, AnnotationInfo &)> f) {for(auto& d: lb.get_lexical_daughters()) f(d, lb.get_annotations());}
  inline void process(function<void(BinaryDaughter &, AnnotationInfo &)> f) { for(auto& d: lb.get_binary_daughters()) f(d, lb.get_annotations()); }

  inline void process(function<void(Best &, LBEdge &, const LexicalDaughter &)> f) {for(const auto& d: lb.get_lexical_daughters()) f(lb.get_best(), lb, d);}
  inline void process(function<void(Best &, LBEdge &, const BinaryDaughter &)> f) {for(const auto& d: lb.get_binary_daughters()) f(lb.get_best(), lb, d);}

  inline void process(function<void(ProbaModel &, const LexicalDaughter &)> f) {for(const auto& d: lb.get_lexical_daughters()) f(lb.get_prob_model(), d);}
  inline void process(function<void(ProbaModel &, const BinaryDaughter &)> f) {for(const auto& d: lb.get_binary_daughters()) f(lb.get_prob_model(), d);}
  inline void process(function<void(ProbaModel &, LexicalDaughter &)> f) {for(auto& d: lb.get_lexical_daughters()) f(lb.get_prob_model(), d);}
  inline void process(function<void(ProbaModel &, BinaryDaughter &)> f) {for(auto& d: lb.get_binary_daughters()) f(lb.get_prob_model(), d);}

  inline void process(function<void(Best &, const LexicalDaughter &)> f) {for(const auto& d: lb.get_lexical_daughters()) f(lb.get_best(), d);}
  inline void process(function<void(Best &, const BinaryDaughter &)> f) {for(const auto& d: lb.get_binary_daughters()) f(lb.get_best(), d);}
  inline void process(function<void(Best &, LexicalDaughter &)> f) {for(auto& d: lb.get_lexical_daughters()) f(lb.get_best(), d);}
  inline void process(function<void(Best &, BinaryDaughter &)> f) {for(auto& d: lb.get_binary_daughters()) f(lb.get_best(), d);}

  void _apply() const {}
  template<typename Function, typename... OtherFunctions>
  void _apply(Function&& f, OtherFunctions&&... o) {process(toFunc(f));_apply(o...);}
  template<typename Function, typename... OtherFunctions>
  void _apply(Function&& f, OtherFunctions&&... o) const {process(toFunc(f));_apply(o...);}

  template<typename... Functions>  void apply(Functions&&... f) const {if (this->is_opened()) _apply(f...);}
  template<typename... Functions>  void apply(Functions&&... f)       {if (this->is_opened()) _apply(f...);}

  template<typename... Functions>  void apply_u(Functions&&... f) const {if (u.is_opened()) _apply(f...);}
  template<typename... Functions>  void apply_u(Functions&&... f)       {if (u.is_opened()) _apply(f...);}

  template<typename... Functions>  void apply_lb(Functions&&... f) const {if (lb.is_opened()) _apply(f...);}
  template<typename... Functions>  void apply_lb(Functions&&... f)       {if (lb.is_opened()) _apply(f...);}

  template<class OPEP>
  friend std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge);
  void dump(std::ostream & out) const;
};

template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const AnnotatedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const BasePackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const UPackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const LBPackedEdge<OPEP>& edge) {edge.dump(out); return out; }
template<class OPEP>
inline std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge) {edge.dump(out); return out; }


#endif /*PACKEDEDGE_H_*/
