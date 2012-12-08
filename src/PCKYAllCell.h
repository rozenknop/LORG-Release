// -*- mode: c++ -*-
#ifndef PCKYALLCELL_H
#define PCKYALLCELL_H

#include "Word.h"

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"
#include "grammars/AnnotatedLabelsInfo.h"


#include <cassert>
#include <cstring>


#include <numeric>
#include <algorithm>
#include <functional>


using std::function;

/**
   \class PCKYAllCell
   \brief represents a cell in the chart for the 2 stage parser
   \note maybe fusion with PCKYBestCell
 */
template<class Types>
class PCKYAllCell {

public:
  typedef typename Types::Edge Edge;
  typedef typename Types::Cell Cell;
  typedef typename Types::BRule BinaryRule;
  typedef typename Types::URule UnaryRule;
  typedef typename Types::LRule LexicalRule;
  typedef typename Types::UnaryDaughter UnaryDaughter;
  typedef typename Types::BinaryDaughter BinaryDaughter;
  typedef typename Types::LexicalDaughter LexicalDaughter;

public:

private:
  /**
     \brief Simple constructor
     \note this cell is not initialised: trying
     to use the created cell will result in segfault !
     You have to call init first
  */
  PCKYAllCell();
  PCKYAllCell(const PCKYAllCell<Types> &);

public:
  /**
     \brief destructor
   */
  ~PCKYAllCell();


  /**
     \brief initialise the cell
     \param cl true if closed
   */
  void init(bool cl, unsigned begin, unsigned end, Edge * first_edge, bool top);
  void reinit(bool cl);


  /**
   * @brief reserve space for binary daughters
   *
   * @param  a vector of int : number of binary daughters to allocate for each edge
   * @return void
   **/
  void reserve_binary_daughters(const std::vector<int> &);
  
  /**
     \brief insert a candidate edge in the cell from application of a binary rule
     \param
  */
  void process_candidate(Edge & left, Edge & right, const BinaryRule*, double LR_inside);


  /**
     \brief insert a candidate edge in the cell from application of a unary rule
     \param
  */
  void process_candidate(const UnaryRule *, double);

  /**
     \brief test if there's an edge of given lhs in the cell
     \param label a lhs label
     \return true if label is in the cell
  */
  bool exists_edge(int label) const;


  /**
     \brief test if there's an edge in the cell
     \return false if there exists at least an edge in the cell
  */
  bool is_empty() const;


/**
     \brief access an edge by its lhs
     \param i the label of the edge
     \return the edge with i as lhs
   */
  Edge * get_edge_ptr(unsigned i);


  /**
     \brief access an edge by its lhs
     \param i the label of the edge
     \return the edge with i as lhs
  */
  const Edge& get_edge(int i) const;
  Edge& get_edge(int i);

  /**
     \brief access
     \return true if the cell is closed
  */
  bool is_closed() const;

  inline bool get_top() const { return top; }
  inline unsigned get_begin() const { return begin; }
  inline unsigned get_end() const { return end; }


  /**
     \brief Output operator
     \param out the ostream to write on
     \param cell the cell object to write
     \return the used ostream
  */
  template<class O>
  friend std::ostream& operator<<(std::ostream& out, const PCKYAllCell<O>& cell);


  /**
    \brief reset probabilities of all the edges in the cell to 0.0
    \note reset to a value passed as a parameter
  */
  void reset_probabilities();

  void compute_inside_probabilities();
  void compute_outside_probabilities();
  void adjust_inside_probability();

  void add_word(const Word & word);

  void beam(const std::vector<double>& priors, double threshold);
  void beam(double threshold);
  void beam(double threshold, double sent_prob);
  void beam_huang(double threshold, double sent_prob);

  static void set_max_size(unsigned size);

  void clean();

  void clear();

  void change_rules_resize(const AnnotatedLabelsInfo& next_annotations, const std::vector<std::vector<std::vector<unsigned> > >& annot_descendants_current);
  void change_rules_resize(unsigned new_size, unsigned fienr_idx);



  /**
   * \brief apply a list of functions on each edge of this cell
   */
  template<typename... Function>
  void apply_on_edges(Function&&... args) {
    for(unsigned i=0; i<get_max_size(); ++i) {/*std::cout << "edge " << i << std::endl ;*/
      if(exists_edge(i)) get_edge(i).apply(args...);
    }
  }

  inline static unsigned get_max_size() { return max_size; }

private:
  //   std::vector<Edge> edges;
  Edge * edges;
  bool closed;
  unsigned begin;
  unsigned end;
  bool top;

  static unsigned max_size;
};

template<class Edge>
Edge * begin(Edge e[]) { return e; }

template<class Edge>
Edge * end(Edge e[]) { return (e+Edge::Cell::get_max_size()); }

template<class Types>
inline
bool PCKYAllCell<Types>::exists_edge(int label) const
{
  assert(label >= 0);
  assert(label < (int) max_size);
  return (not edges[label].is_closed());
}


template<class Types>
inline
bool PCKYAllCell<Types>::is_empty() const
{
  if(closed)
    return true;

  for (unsigned i = 0; i < max_size; ++i)
    {
      if(not edges[i].is_closed())
        return false;
    }
  return true;
}


template<class Types>
inline
void PCKYAllCell<Types>::init(bool cl, unsigned b, unsigned e, typename Types::Edge * first_edge, bool t)
{
  begin = b; end = e; top = t; closed = cl; edges = first_edge;
}

template<class Types>
inline
void PCKYAllCell<Types>::reinit(bool cl)
{
  if(!(closed = cl)) {
    clear();
  }
}

template<class Types>
inline
bool PCKYAllCell<Types>::is_closed() const
{ return closed; }

template<class Types>
inline
typename Types::Edge * PCKYAllCell<Types>::get_edge_ptr(unsigned i)
{
  assert(!closed);
  return edges[i];
}

template<class Types>
inline
const typename Types::Edge& PCKYAllCell<Types>::get_edge(int i) const
{
  //assert(i>=0);
  //assert( i < (int) max_size);
  assert(i>=0 && i < (int) max_size);

  return edges[i];
}

template<class Types>
inline
typename Types::Edge& PCKYAllCell<Types>::get_edge(int i)
{
  //assert(i>=0);
  //assert( i < (int) max_size);
  assert(i>=0 && i < (int) max_size);

  return edges[i];
}


template<class Types>
inline
void PCKYAllCell<Types>::set_max_size(unsigned size)
{
  max_size =  size;
}
#endif //PCKYALLCELL_H
