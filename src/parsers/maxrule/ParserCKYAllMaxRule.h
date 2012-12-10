// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVAR_H_
#define _PARSERCKYALLMAXVAR_H_

//#include "ParserCKY.h"
#include "parsers/ParserCKYAll.h"
#include "utils/lorg_functional.h"

template<class Types>
class ParserCKYAllMaxRule : public ParserCKYAll_Impl<Types>
{
public:
  typedef typename ParserCKYAll_Impl<Types>::AGrammar AGrammar;
  typedef typename Types::Chart Chart;
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::EdgeProbability ProbaModel;

    ParserCKYAllMaxRule(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& priors,
                        double beam_threshold,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam_length, int stubborn)
    : ParserCKYAll_Impl<Types>(cgs,
                               priors,
                               beam_threshold,
                               annot_descendants_,
                               accurate_, min_beam_length, stubborn)
    {}



   /**
   * @brief calculate the chart specific rule probabilities for all packed edges in all cells
   *
   * @return void
   **/
   void calculate_maxrule_probabilities()
   {
     this->chart->opencells_apply_bottom_up(
      [](Cell & cell)
      {
//       std::cout << "filling cell " << &cell << " : ======================================================" << cell << std::endl;
        cell.apply_on_edges (toFunc(&ProbaModel::update_lexical),
                             toFunc (&ProbaModel::update_binary));
        cell.apply_on_edges (toFunc(&ProbaModel::update_unary),
                             toFunc (&ProbaModel::finalize));
//       std::cout << "best filled for cell " << &cell << " : " << cell << std::endl;
      }
    );
  }
};

#endif /* _PARSERCKYALLMAXVAR_H_ */
