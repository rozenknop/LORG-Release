// -*- mode: c++ -*-
#ifndef PACKEDEDGEDAUGHTERS_H
#define PACKEDEDGEDAUGHTERS_H

#include "utils/tick_count.h"

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"
#include "AnnotationInfo.h"

class Word;
//class PCKYAllCell;

/**
   \class PackedEdgeDaughters
   \brief represents a branching of possible daughters
*/

class PackedEdgeDaughters
{
protected:
  const AnnotatedRule * rule;

public:
  PackedEdgeDaughters(const AnnotatedRule * r = 0) : rule(r) {};
  ~PackedEdgeDaughters() {};

  inline bool is_unary() const {return get_rule()->is_unary();}
  inline bool is_binary() const {return get_rule()->is_binary();}
  inline bool is_lexical() const {return get_rule()->is_lexical();}

// protected:
  inline const AnnotatedRule * get_rule() const {return rule;}
  inline void set_rule(const AnnotatedRule * r) {rule = r;}

};



template<class Rule>
class TypedRulePackedEdgeDaughters : public PackedEdgeDaughters
{
public:
  TypedRulePackedEdgeDaughters (const Rule * r) : PackedEdgeDaughters(r) {}
  inline const Rule * get_rule() const {return (Rule*) rule;}
};


/**
   \class BinaryPackedEdgeDaughters
   \brief represents a binary branching of possible daughters + a binary rule
*/
template<class Types>
class BinaryPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::BRule>
{
public:
  typedef TypedRulePackedEdgeDaughters<typename Types::BRule> Parent;
  typedef typename Types::Cell Cell;
  typedef typename Types::AEdge AEdge;
  typedef typename Types::PEdge PEdge;
  typedef typename Types::LBEdge LBEdge;
  typedef typename Types::UEdge UEdge;
  typedef typename Types::AEdge * EdgePtr ;
  typedef typename Types::BRule Rule;

protected:
  EdgePtr left;
  EdgePtr right;
//   friend class Types::Edge ;
public:
  BinaryPackedEdgeDaughters(AEdge& le, AEdge& ri, const typename Types::BRule * ru) :
    Parent(ru), left(&le),right(&ri)
  {};

  ~BinaryPackedEdgeDaughters() {};

//   inline const Edge& left_daughter() const  {return *left;}
//   inline const Edge& right_daughter() const {return *right;}
  inline AEdge& left_daughter() const {return *left;}
  inline AEdge& right_daughter() const {return *right;}
  inline PEdge& left_pdaughter() const {return (PEdge &) *left;}
  inline PEdge& right_pdaughter() const {return (PEdge &) *right;}

  inline bool operator==(const BinaryPackedEdgeDaughters& other)
  {
    return Parent::rule == other.rule && left == other.left && right ==other.right;
  }
  inline bool points_towards_invalid_edges() const
  {
    return left->is_closed() or right->is_closed() ;
  }
  
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(Parent::rule != NULL);
    Parent::get_rule()->update_inside_annotations(annotations.inside_probabilities.array,
                                                  left-> get_annotations().inside_probabilities.array,
                                                  right->get_annotations().inside_probabilities.array);
  }
  
  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    Parent::get_rule()->update_outside_annotations(annotations.outside_probabilities.array,
                                        left->get_annotations().inside_probabilities.array,
                                        right->get_annotations().inside_probabilities.array,
                                        left->get_annotations().outside_probabilities.array,
                                        right->get_annotations().outside_probabilities.array);
  }
};


/**
   \class UnaryPackedEdgeDaughters
   \brief represents a unary branching (!) + a unary rule
*/
template<class Types>
class UnaryPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::URule>
{
public:
  typedef TypedRulePackedEdgeDaughters<typename Types::URule> Parent;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::AEdge AEdge;
  typedef typename Types::LBEdge LBEdge;
  typedef typename Types::UEdge UEdge;
  typedef typename Types::AEdge * EdgePtr ;
  typedef typename Types::URule Rule;

protected:
  EdgePtr dtr;

public:
//   inline UnaryPackedEdgeDaughters & operator=(UnaryPackedEdgeDaughters<Types> && o) { *this = std::move(o); return *this; }
  
  UnaryPackedEdgeDaughters(AEdge & le, const Rule * ru) :
  Parent(ru), dtr(&le)
  {};

  ~UnaryPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return false;}
  //   inline const Edge& left_daughter() const  {return *left;}
  inline AEdge& daughter() const {return *dtr;}
  inline LBEdge& lbdaughter() const {return (LBEdge &)(*dtr);}

  inline bool points_towards_invalid_edges() const
  {
    return dtr->is_closed();
  }
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(Parent::rule != NULL);
    Parent::get_rule()->update_inside_annotations(annotations.inside_probabilities.array,
                                                  dtr->get_annotations().inside_probabilities.array);
  }
  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    Parent::get_rule()->update_outside_annotations(annotations.outside_probabilities.array,
                                                   dtr->get_annotations().outside_probabilities.array);
  }
};


/**
   \class LexicalPackedEdgeDaughters
   \brief represents a zero-ary branching (!) + a lexical rule
*/
template<class Types>
class LexicalPackedEdgeDaughters : public TypedRulePackedEdgeDaughters<typename Types::LRule>
{
public:
  typedef typename Types::LRule Rule;
  typedef TypedRulePackedEdgeDaughters<typename Types::LRule> Parent;

  protected:
  const Word* word;

  double relaxation;

public:
  LexicalPackedEdgeDaughters(const Rule * ru, const Word* w) :
  Parent(ru), word(w),
    relaxation(0)
  {};

  ~LexicalPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return true;}
  inline const Word* get_word() const {return word;}
  
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
//     BLOCKTIMING("LexicalPackedEdgeDaughters - update_inside_annotations");
    assert(Parent::get_rule() != NULL);
    Parent::get_rule()->update_inside_annotations(annotations.inside_probabilities.array);
  }
};


#endif //PACKEDEDGEDAUGHTERS
