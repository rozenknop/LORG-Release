// -*- mode: c++ -*-
#ifndef CHARTCKY_H_
#define CHARTCKY_H_

#include "utils/PtbPsTree.h"
#include "utils/LorgConstants.h"
#include "utils/SymbolTable.h"

#include <vector>
#include <Bracketing.h>

#ifdef USE_THREADS
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#endif

#include <functional>
#include <iostream>
using std::ostream;
using std::endl;

/**
  \class ChartCKY
  \brief represents a chart of cells
*/
template<class TCell, class MyWord>
class ChartCKY
{
public:
    typedef TCell Cell;

private:
  TCell ** chart; ///< the chart itself
  unsigned size;     ///< the size of the chart
  const std::vector< MyWord >& sentence;
  const std::vector<bracketing>& brackets;
  #ifdef USE_THREADS
  std::vector<Cell *> vcells;
  #endif

  // prevents unwanted conversions
  ChartCKY(const ChartCKY&);
  ChartCKY& operator=(const ChartCKY&);

public:
  ~ChartCKY();

  /**
     \brief constructor with initialisation
     \param sentence the sentence to create the chart
     \param grammar_size  the number non-terminals in the grammar
     \param brackets chunks
  */
  ChartCKY(const std::vector< MyWord >& sentence, unsigned grammar_size, const std::vector<bracketing>& brackets);

  /**
     \brief get the size of the chart
     \return the size of the chart
  */
  unsigned get_size() const;


  /**
     \brief access a cell of the chart by its coordinates
     \param start starting point of the cell's span
     \param end   ending point of the cell's span
     \return a cell (may segfault if coordinates are out of bounds)
  */

  TCell& access(unsigned start, unsigned end) const;

  TCell& get_root() const;

  PtbPsTree* get_best_tree(int start_symbol, unsigned k, bool always_output_forms, bool output_annotations) const;

  double get_score(int start_symbol, unsigned k) const;


  void init(const std::vector< MyWord >& sentence);

  void reset_probabilities();

  bool has_solution(int symb, unsigned i) const;

  void clear();

  void prepare_retry();

  bool is_valid(int start_symbol) const;

  ostream & to_stream(ostream & s) const;

  void opencells_apply( std::function<void(Cell &)> f );
  void opencells_apply_nothread( std::function<void(Cell &)> f );
  void opencells_apply_bottom_up( std::function<void(Cell &)> f, unsigned min_span=0 );
  void opencells_apply_bottom_up_nothread( std::function<void(Cell &)> f, unsigned min_span=0 );
  void opencells_apply_top_down( std::function<void(Cell &)> f );
  void opencells_apply_top_down_nothread( std::function<void(Cell &)> f );

  std::ostream & operator>>(std::ostream & out) { opencells_apply_bottom_up([out](TCell & cell){return out << cell << endl; }); return out; }
};


#endif /*CHARTCKY_H_*/
