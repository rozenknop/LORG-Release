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
template<class Types>
class ChartCKY
{
public:
  typedef typename Types::Cell Cell;
  typedef typename Types::Edge Edge;
  typedef typename Types::ChartWord MyWord;

private:
  Cell * the_cells; ///< the chart itself
  Edge * the_edges; ///< the edges of the chart
//   Cell ** chart; ///< pointers on each column of the chart
  unsigned size;     ///< the size of the chart (width)
  unsigned nb_cells; ///< number of cells in the chart
  const std::vector< MyWord >& sentence;
  const std::vector<bracketing>& brackets;

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
  inline unsigned get_size() const;


  /**
     \brief access a cell of the chart by its coordinates
     \param start starting point of the cell's span
     \param end   ending point of the cell's span
     \return a cell (may segfault if coordinates are out of bounds)
  */

  inline const Cell& access(unsigned start, unsigned end) const;
  inline Cell& access(unsigned start, unsigned end);

  inline const Cell& get_root() const;
  inline Cell& get_root();

  inline PtbPsTree* get_best_tree(int start_symbol, unsigned k, bool always_output_forms, bool output_annotations, bool unary_start=true) const;

  inline double get_score(int start_symbol, unsigned k, bool unary_start=true) const;


  inline void reset_probabilities();

  inline bool has_solution(int symb, unsigned i, bool unary_start=true) const;

  inline void clear();

  inline void prepare_retry();

  inline bool is_valid(int start_symbol) const;

  void dump(ostream & s) const ;

  inline void opencells_apply( std::function<void(Cell &)> f );
  inline void opencells_apply_nothread( std::function<void(Cell &)> f );
  inline void opencells_apply_bottom_up( std::function<void(Cell &)> f, unsigned min_span=0 );
  inline void opencells_apply_bottom_up_nothread( std::function<void(Cell &)> f, unsigned min_span=0 );
  inline void opencells_apply_top_down( std::function<void(Cell &)> f );
  inline void opencells_apply_top_down_nothread( std::function<void(Cell &)> f );
  inline void opencells_apply_top_down_nothread( std::function<void(const Cell &)> f ) const;
  inline void opencells_apply_left_right( std::function<void(const Cell &)> before_left,
                                          std::function<void(const Cell &)> before_right,
                                          std::function<void(const Cell &)> after_daughters);
  inline void opencells_apply_left_right_rec( std::function<void(const Cell &)> before_left,
                                              std::function<void(const Cell &)> before_right,
                                              std::function<void(const Cell &)> after_daughters,
                                              unsigned span, unsigned beg, bool go_left
                                            );
};

template<class OPEP>
inline std::ostream & operator<<(std::ostream & out, const ChartCKY<OPEP> & chart) { chart.dump(out); return out; }

template <class Types>
inline bool ChartCKY<Types>::is_valid(int start_symbol) const
{
  return !get_root().is_closed() && get_root().exists_uedge(start_symbol);
}

template<class Types>
inline typename Types::Cell& ChartCKY<Types>::get_root()
{
  return access(0,size-1);
}

template<class Types>
inline const typename Types::Cell& ChartCKY<Types>::get_root() const
{
  return access(0,size-1);
}

template<class Types>
inline
unsigned ChartCKY<Types>::get_size() const
{
  return size;
}

template<class Types>
inline
typename Types::Cell& ChartCKY<Types>::access(unsigned start, unsigned end)
{
  //   assert(start <= end);
  //   assert(end < size);
  //     return chart[start][end-start];
  //     std::cout << "access("<<start<<","<<end<<") = "<<start + ( (end-start)*(2*size - (end-start) + 1))/2 << std::endl; std::cout.flush();
  return the_cells[start + ( (end-start)*(2*size - (end-start) + 1))/2 ];
}

template<class Types>
inline
const typename Types::Cell& ChartCKY<Types>::access(unsigned start, unsigned end) const
{
  //   assert(start <= end);
  //   assert(end < size);
  //     return chart[start][end-start];
  return the_cells[start + ( (end-start)*(2*size - (end-start) + 1))/2 ];
}


#endif /*CHARTCKY_H_*/
