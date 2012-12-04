// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVARKB_H_
#define _PARSERCKYALLMAXVARKB_H_

#include "edges/MaxRuleProbabilityKB.h"
#include "ParserCKYAllMaxVar.h"


class ParserCKYAllMaxRuleKB : public ParserCKYAllMaxRule<MaxRuleKBTypes>
{
 private:
  unsigned k;

 public:
  ParserCKYAllMaxRuleKB(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& p, double b_t,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam, int stubborn, unsigned k_);

  ~ParserCKYAllMaxRuleKB() {};

  void extract_solution();

 private:
  void initialise_candidates();

  void extend_all_derivations();
};

#endif /* _PARSERCKYALLMAXVARKB_H_ */
