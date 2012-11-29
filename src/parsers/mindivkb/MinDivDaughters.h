#ifndef _MINDIVDAUGHTER_H_
#define _MINDIVDAUGHTER_H_

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"

#include "edges/PackedEdgeDaughters.h"
#include "MinDivTypes.h"



struct MinDivBRule : public BRuleC2f {
  double update_outside_annotations_return_marginal(const std::vector<double>& up_out,
                                                    const std::vector<double>& left_in,
                                                    const std::vector<double>& right_in,
                                                    std::vector<double>& left_out,
                                                    std::vector<double>& right_out) const;
};

struct MinDivURule : public URuleC2f {
  double update_outside_annotations_return_marginal(const std::vector<double>& up,
                                                    const std::vector<double>& in_left,
                                                    std::vector<double>& out_left) const;
};

struct MinDivLRule : public LexicalRuleC2f {
  double update_outside_annotations_return_marginal(const std::vector<double>& up) const;
};






struct MinDivEdgeDaughterProbability
{
  /// marginal probability on p (annotated forest)
  double mp;
  /// rule probability on q
  double q;
};


class MinDivBinaryDaughter : public BinaryPackedEdgeDaughters<MinDivKBTypes>, public MinDivEdgeDaughterProbability
{
public:
  MinDivBinaryDaughter(Cell *le, Cell *ri, const MinDivBRule * ru) : BinaryPackedEdgeDaughters<MinDivKBTypes>(le,ri,ru) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) {BinaryPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}
  inline void update_outside_annotations(AnnotationInfo & annotations) {BinaryPackedEdgeDaughters<MinDivKBTypes>::update_outside_annotations(annotations);}

  inline void outside_and_marginal(AnnotationInfo & annotations)
  {
    auto & leftannot = left_daughter()->get_edge(get_rule()->get_rhs0()).get_annotations();    
    auto & rightannot= right_daughter()->get_edge(get_rule()->get_rhs1()).get_annotations();
    mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                                leftannot.inside_probabilities.array,
                                                                rightannot.inside_probabilities.array,
                                                                leftannot.outside_probabilities.array,
                                                                rightannot.outside_probabilities.array)
    / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
  }
};

class MinDivUnaryDaughter : public UnaryPackedEdgeDaughters<MinDivKBTypes>, public MinDivEdgeDaughterProbability
{
public:
  MinDivUnaryDaughter(Cell *le, const MinDivURule * ru) : UnaryPackedEdgeDaughters<MinDivKBTypes>(le,ru) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) {UnaryPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}
  inline void update_outside_annotations(AnnotationInfo & annotations) {UnaryPackedEdgeDaughters<MinDivKBTypes>::update_outside_annotations(annotations);}
  
  inline void outside_and_marginal(AnnotationInfo & annotations)
  {
    auto & leftannot = left_daughter()->get_edge(get_rule()->get_rhs0()).get_annotations();
    mp = RH::rule->update_outside_annotations_return_marginal(annotations.outside_probabilities.array,
                                                              leftannot.inside_probabilities.array,
                                                              leftannot.outside_probabilities_unary_temp.array)
         / MinDivProbabilityKB::get_normalisation_factor();
  }
};

class MinDivLexicalDaughter : public LexicalPackedEdgeDaughters<MinDivKBTypes>, public MinDivEdgeDaughterProbability
{
public:
  MinDivLexicalDaughter(const MinDivLRule * ru, const Word* w) : LexicalPackedEdgeDaughters<MinDivKBTypes>(ru,w) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) {LexicalPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}
  
  inline void outside_and_marginal(AnnotationInfo & annotations)
  {
    mp = get_rule()->update_outside_annotations_return_marginal(annotations.outside_probabilities.array)
    / MinDivProbabilityKB::get_normalisation_factor();
//     std::cout << mp << " de " << MinDivProbabilityKB::get_normalisation_factor() << std::endl;
  }
  
  inline void update_inside(MinDivKBTypes::Edge & edge) {
    edge.get_prob_model().get_inside_prob() += q ;
  };
  
  inline void update_outside(MinDivKBTypes::Edge &) {}    
};

#endif
