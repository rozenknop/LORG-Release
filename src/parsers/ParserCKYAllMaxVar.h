// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMAXVAR_H_
#define _PARSERCKYALLMAXVAR_H_

#include "ParserCKY.h"
#include "ParserCKYAll.h"
#include "utils/lorg_functional.h"

template<class TCell>
class ParserCKYAllMaxRule : public ParserCKYAll_Impl<TCell>
{
public:
    typedef typename ParserCKYAll_Impl<TCell>::AGrammar AGrammar;
    typedef typename ParserCKYAll_Impl<TCell>::Chart Chart;
    typedef typename Chart::Cell Cell;
    typedef typename Cell::Edge Edge;
    typedef typename Edge::ProbaModel ProbaModel;

    ParserCKYAllMaxRule(std::vector<AGrammar*>& cgs,
                        const std::vector<double>& priors,
                        double beam_threshold,
                        const annot_descendants_type& annot_descendants_,
                        bool accurate_, unsigned min_beam_length, int stubborn)
    : ParserCKYAll_Impl<TCell>(cgs,
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
     this->chart->opencells_apply_bottom_up (
      [](Cell & cell)
      {
        cell.apply_on_edges (toFunc(&ProbaModel::update_lexical),
                             toFunc (&ProbaModel::update_binary));
        cell.apply_on_edges (toFunc(&ProbaModel::update_unary),
                             toFunc (&ProbaModel::finalize));
      },
      this->num_cell_threads
    );
  }

};

#endif /* _PARSERCKYALLMAXVAR_H_ */
