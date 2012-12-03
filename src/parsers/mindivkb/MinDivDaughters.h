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
//   MinDivBinaryDaughter(const MinDivBinaryDaughter&&other) : BinaryPackedEdgeDaughters<MinDivKBTypes>(other) {left=other.left; right=other.right; rule=other.rule;}
  inline MinDivBinaryDaughter(Cell *le, Cell *ri, const MinDivKBTypes::BRule * ru) : BinaryPackedEdgeDaughters<MinDivKBTypes>(le,ri,ru) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) const {BinaryPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}
  inline void update_outside_annotations(AnnotationInfo & annotations) const {BinaryPackedEdgeDaughters<MinDivKBTypes>::update_outside_annotations(annotations);}

  inline void outside_and_marginal(AnnotationInfo & annotations);
  inline double tree_log_proba(unsigned left_idx =0, unsigned right_idx = 0) const;
};

class MinDivUnaryDaughter : public UnaryPackedEdgeDaughters<MinDivKBTypes>, public MinDivEdgeDaughterProbability
{
public:
  inline MinDivUnaryDaughter(Cell *le, const MinDivKBTypes::URule * ru) : UnaryPackedEdgeDaughters<MinDivKBTypes>(le,ru) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) const {UnaryPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}
  inline void update_outside_annotations(AnnotationInfo & annotations) const {UnaryPackedEdgeDaughters<MinDivKBTypes>::update_outside_annotations(annotations);}

  inline void outside_and_marginal(AnnotationInfo & annotations);
  inline double tree_log_proba(unsigned left_idx =0) const;
};

class MinDivLexicalDaughter : public LexicalPackedEdgeDaughters<MinDivKBTypes>, public MinDivEdgeDaughterProbability
{
public:
  inline MinDivLexicalDaughter(const MinDivKBTypes::LRule * ru, const Word* w) : LexicalPackedEdgeDaughters<MinDivKBTypes>(ru,w) {}
  inline void update_inside_annotations(AnnotationInfo & annotations) const {LexicalPackedEdgeDaughters<MinDivKBTypes>::update_inside_annotations(annotations);}

  inline void outside_and_marginal(AnnotationInfo & annotations);
  inline double tree_log_proba() const;
};

#endif
