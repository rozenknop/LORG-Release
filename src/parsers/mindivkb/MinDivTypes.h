#ifndef _MINDIVTYPES_H_
#define _MINDIVTYPES_H_

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
  typedef MinDivBRule BRule;
  typedef MinDivURule URule;
  typedef MinDivLRule LRule;
  typedef MinDivBinaryDaughter BinaryDaughter;
  typedef MinDivUnaryDaughter  UnaryDaughter;
  typedef MinDivLexicalDaughter LexicalDaughter;
};

#endif
