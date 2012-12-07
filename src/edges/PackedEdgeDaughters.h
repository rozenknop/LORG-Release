// -*- mode: c++ -*-
#ifndef PACKEDEDGEDAUGHTERS_H
#define PACKEDEDGEDAUGHTERS_H

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"
#include "AnnotationInfo.h"

class Word;
//class PCKYAllCell;

/**
   \class RuleHolder
   \brief a class with an attribute rule
 */
template<class MyRule>
class RuleHolder
{
public:
  typedef MyRule ruletype;
protected:
  const MyRule * rule;
public:
  typedef MyRule rule_type;
  RuleHolder(const MyRule * r = 0) : rule(r) {};
  ~RuleHolder() {};
  inline const MyRule * get_rule() const {return rule;}
  inline void set_rule(const MyRule * r) {rule = r;}
};


/**
   \class PackedEdgeDaughters
   \brief represents a branching of possible daughters
*/

class PackedEdgeDaughters
{

public:
  PackedEdgeDaughters()  {};

  virtual bool is_binary() const =0;
  virtual bool is_lexical() const =0;
  ~PackedEdgeDaughters() {};
};


/**
   \class BinaryPackedEdgeDaughters
   \brief represents a binary branching of possible daughters + a binary rule
*/
template<class Types>
class BinaryPackedEdgeDaughters : public PackedEdgeDaughters, public RuleHolder<typename Types::BRule>
{
public:
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::Edge * EdgePtr ;
  typedef typename Types::BRule Rule;
  typedef RuleHolder<Rule> RH;
protected:
  EdgePtr left;
  EdgePtr right;
//   friend class Types::Edge ;
public:
//   inline BinaryPackedEdgeDaughters & operator=(BinaryPackedEdgeDaughters<Types> && o) { *this = std::move(o); return *this; }
  BinaryPackedEdgeDaughters(Edge& le, Edge& ri, const typename Types::BRule * ru) :
    PackedEdgeDaughters(), RH(ru), left(&le),right(&ri)
  {};
//   BinaryPackedEdgeDaughters(BinaryPackedEdgeDaughters&& o) : PackedEdgeDaughters(), RH(o), left(o.left),right(o.right) {}
//   BinaryPackedEdgeDaughters(const BinaryPackedEdgeDaughters& o) : PackedEdgeDaughters(), RH(o), left(o.left),right(o.right) {}

  ~BinaryPackedEdgeDaughters() {};

  inline bool is_binary() const {return true;}
  inline bool is_lexical() const {return false;}

//   inline const Edge& left_daughter() const  {return *left;}
//   inline const Edge& right_daughter() const {return *right;}
  inline Edge& left_daughter() const {return *left;}
  inline Edge& right_daughter() const {return *right;}

  inline bool operator==(const BinaryPackedEdgeDaughters& other)
  {
    return RH::rule == other.rule && left == other.left && right ==other.right;
  }
  inline bool points_towards_invalid_edges() const
  {
    return left->is_closed() or right->is_closed() ;
  }
  
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(RH::rule != NULL);
    RH::rule->update_inside_annotations(annotations.inside_probabilities.array,
                                        left->get_annotations().inside_probabilities.array,
                                        right->get_annotations().inside_probabilities.array);
  }
  
  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    RH::rule->update_outside_annotations(annotations.outside_probabilities.array,
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
class UnaryPackedEdgeDaughters : public PackedEdgeDaughters, public RuleHolder<typename Types::URule>
{
public:
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::Edge * EdgePtr ;
  typedef typename Types::URule Rule;
  typedef RuleHolder<Rule> RH;
protected:
  EdgePtr left;

public:
//   inline UnaryPackedEdgeDaughters & operator=(UnaryPackedEdgeDaughters<Types> && o) { *this = std::move(o); return *this; }
  
  UnaryPackedEdgeDaughters(Edge & le, const Rule * ru) :
    PackedEdgeDaughters(), RH(ru), left(&le)
  {};

  ~UnaryPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return false;}
//   inline const Edge& left_daughter() const  {return *left;}
  inline Edge& left_daughter() const {return *left;}
  
  inline bool points_towards_invalid_edges() const
  {
    return left->is_closed();
  }
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(RH::rule != NULL);
    RH::rule->update_inside_annotations(annotations.inside_probabilities_unary_temp.array,
                                        left->get_annotations().inside_probabilities.array);
  }
  inline void update_outside_annotations(AnnotationInfo & annotations) const
  {
    RH::rule->update_outside_annotations(annotations.outside_probabilities.array,
                                         left->get_annotations().outside_probabilities_unary_temp.array);
  }
};


/**
   \class LexicalPackedEdgeDaughters
   \brief represents a zero-ary branching (!) + a lexical rule
*/
template<class Types>
class LexicalPackedEdgeDaughters : public PackedEdgeDaughters, public RuleHolder<typename Types::LRule>
{
public:
  typedef typename Types::LRule Rule;
  typedef RuleHolder<Rule> RH;

  protected:
  const Word* word;

  double relaxation;

public:
  LexicalPackedEdgeDaughters(const Rule * ru, const Word* w) :
  PackedEdgeDaughters(), RH(ru), word(w),
    relaxation(0)
  {};

  ~LexicalPackedEdgeDaughters() {};

  inline bool is_binary() const {return false;}
  inline bool is_lexical() const {return true;}
  inline const Word* get_word() const {return word;}
  
  inline void update_inside_annotations(AnnotationInfo & annotations) const {
    assert(RH::rule != NULL);
    RH::rule->update_inside_annotations(annotations.inside_probabilities.array);
  }
};


#endif //PACKEDEDGEDAUGHTERS
