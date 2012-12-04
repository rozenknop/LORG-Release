#ifndef _MINDIVTYPES_H_
#define _MINDIVTYPES_H_
#include <rules/BRuleC2f.h>
#include <rules/LexicalRuleC2f.h>
#include <rules/URuleC2f.h>

class MinDivProbabilityKB;
class MinDivEdgeDaughterProbability;
class MinDivBRule;
class MinDivURule;
class MinDivLRule;
class MinDivBinaryDaughter;
class MinDivUnaryDaughter;
class MinDivLexicalDaughter;


struct MinDivKBTypes {
  typedef MinDivProbabilityKB EdgeProbability ;
  typedef MinDivEdgeDaughterProbability EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef PackedEdge< MinDivKBTypes > Edge ;
  typedef PCKYAllCell< MinDivKBTypes > Cell ;
  typedef ChartCKY< MinDivKBTypes > Chart ;
//   typedef MinDivBRule BRule;
//   typedef MinDivURule URule;
//   typedef MinDivLRule LRule;
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef MinDivBinaryDaughter BinaryDaughter;
  typedef MinDivUnaryDaughter  UnaryDaughter;
  typedef MinDivLexicalDaughter LexicalDaughter;
};

#endif
