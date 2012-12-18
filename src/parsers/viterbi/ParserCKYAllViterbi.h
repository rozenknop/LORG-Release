// -*- mode: c++ -*-
#ifndef _PARSERCKYALLVITERBI_H_
#define _PARSERCKYALLVITERBI_H_

#include "parsers/ParserCKYAll.h"
#include "ViterbiProbability.h"
#include "emptystruct.h"



class ParserCKYAllViterbi : public ParserCKYAll_Impl<ViterbiTypes>
{
public:
  ParserCKYAllViterbi(std::vector<AGrammar*>& cgs,
                      const std::vector<double>& p, double b_t,
                      const std::vector< std::vector<std::vector< std::vector<unsigned> > > >& annot_descendants_, bool accurate_, unsigned min_beam, int stubborn);

  virtual ~ParserCKYAllViterbi() { delete fine_grammar; fine_grammar = NULL;};

  void extract_solution();
  const AGrammar& get_fine_grammar() const;

private: // attributes
  AGrammar* fine_grammar; ///< the grammar to be used to extract the solution
};


#endif /* _PARSERCKYALLVITERBI_H_ */
