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
template<class MyEdge>
class PCKYAllCell {
public:
  typedef  MyEdge Edge;
public:
  /**
     \brief Simple constructor
     \note this cell is not initialised: trying
     to use the created cell will result in segfault !
     You have to call init first
  */
  PCKYAllCell() : edges(nullptr), closed(true) {};

  /**
     \brief Constructor
     \param cl true if closed
  */
  PCKYAllCell(bool cl);

  /**
     \brief destructor
   */
  ~PCKYAllCell();


  /**
     \brief initialise the cell
     \param cl true if closed
   */
  void init(bool cl, unsigned begin, unsigned end, bool top);
  void reinit(bool cl);
  

  /**
     \brief insert a candidate edge in the cell from application of a binary rule
     \param
  */
  void process_candidate(PCKYAllCell* left, PCKYAllCell* right, const BRuleC2f*, double LR_inside);


  /**
     \brief insert a candidate edge in the cell from application of a unary rule
     \param
  */
  void process_candidate(const URuleC2f *, double);

  /**
     \brief test if there's an edge of given lhs in the cell
     \param label a lhs label
     \return true if label is in the cell
  */
  bool exists_edge(int label) const;


  /**
     \brief test if there's an edge in the cell
     \return fasle if there exists at least an edge in the cell
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

  /**
     \brief compute the best viterbi derivations for a cell
  */
  void compute_best_viterbi_derivation(const AnnotatedLabelsInfo& symbol_map);

  void add_word(const Word & word);

  void beam(const std::vector<double>& priors, double threshold);
  void beam(double threshold);
  void beam(double threshold, double sent_prob);
  void beam_huang(double threshold, double sent_prob);

  static void set_max_size(unsigned size);

  void clean();
  void clean_binary_daughters();

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

  unsigned get_max_size() const { return max_size; }
  
private:
  Edge ** edges;
  bool closed;
  unsigned begin;
  unsigned end;
  bool top;

  static unsigned max_size;
};


template<class PEProbability>
inline
bool PCKYAllCell<PEProbability>::exists_edge(int label) const
{
  assert(label >= 0);
  assert(label < (int) max_size);
  return (edges[label] != nullptr);
}


template<class PEProbability>
inline
bool PCKYAllCell<PEProbability>::is_empty() const
{
  if(closed)
    return true;

  for (unsigned i = 0; i < max_size; ++i)
    {
      if(edges[i])
        return false;
    }
  return true;
}




template<class PEProbability>
inline
void PCKYAllCell<PEProbability>::init(bool cl, unsigned b, unsigned e, bool t)
{
  begin = b; end = e; top = t;
  reinit(cl);
}

template<class PEProbability>
inline
void PCKYAllCell<PEProbability>::reinit(bool cl)
{
  if(!(closed = cl)) {
    edges =  new Edge * [max_size];
    memset(edges, 0, max_size * sizeof(Edge*));
    //   for(unsigned i = 0; i < max_size;++i)
    //     edges[i]=NULL;
  }
}

template<class MyEdge>
inline
bool PCKYAllCell<MyEdge>::is_closed() const
{ return closed; }

template<class MyEdge>
inline
MyEdge * PCKYAllCell<MyEdge>::get_edge_ptr(unsigned i)
{
  assert(!closed);
  return edges[i];
}

template<class MyEdge>
inline
const MyEdge& PCKYAllCell<MyEdge>::get_edge(int i) const
{
  //assert(i>=0);
  //assert( i < (int) max_size);
  assert(i>=0 && i < (int) max_size);

  return *edges[i];
}

template<class MyEdge>
inline
MyEdge& PCKYAllCell<MyEdge>::get_edge(int i)
{
  //assert(i>=0);
  //assert( i < (int) max_size);
  assert(i>=0 && i < (int) max_size);

  return *edges[i];
}

template<class MyEdge>
inline
void PCKYAllCell<MyEdge>::add_word(const Word & word)
{
  typedef typename MyEdge::LexicalDaughters LDaughters;

  for(std::vector<const MetaProduction*>::const_iterator it(word.get_rules().begin());
      it != word.get_rules().end(); ++it) {
    // std::cout <<*(static_cast<const LexicalRule*>(*it)) << std::endl;

    int tag = (*it)->get_lhs();

    Edge ** e = &edges[tag];

    if(*e) {
      (*e)->add_daughters(static_cast<const LexicalRuleC2f*>(*it), &word);
    }
    else {
      *e = new Edge(LDaughters(static_cast<const LexicalRuleC2f*>(*it), &word));
    }

    (*e)->get_annotations().inside_probabilities.array[0] += static_cast<const LexicalRuleC2f*>(*it)->get_probability()[0];
  }
}



template<class MyEdge>
inline
void PCKYAllCell<MyEdge>::set_max_size(unsigned size)
{
  max_size =  size;
}

template<class MyEdge>
unsigned PCKYAllCell<MyEdge>::max_size = 0;


template<class MyEdge>
PCKYAllCell<MyEdge>::~PCKYAllCell()
{
  if(!closed) {
    for(unsigned i = 0; i < max_size;++i) {
      delete edges[i];
      //      edges[i] = NULL;
    }
  }
  delete[] edges;
  //  edges = NULL;
}

template<class MyEdge>
void PCKYAllCell<MyEdge>::process_candidate(PCKYAllCell<MyEdge>* left,
                                            PCKYAllCell<MyEdge>* right,
                                            const BRuleC2f* rule,
                                            double LR_inside)
{
  MyEdge ** e = &edges[rule->get_lhs()];

  if(*e)
    (*e)->add_daughters(left,right,rule);
  else {
    *e = new MyEdge(typename MyEdge::BinaryDaughters(left,right,rule));
  }


  (*e)->get_annotations().inside_probabilities.array[0] += LR_inside * rule->get_probability()[0][0][0];
}

template<class MyEdge>
void PCKYAllCell<MyEdge>::process_candidate(const URuleC2f* rule, double L_inside)
{
  assert(rule);
  assert(rule->get_probability().size() > 0);


  MyEdge ** e = &edges[rule->get_lhs()];

  if(*e)  {
    (*e)->add_daughters(this,rule);
  }
  else {
    //std::cout <<" adding a new edge " << *rule << std::endl;
    *e = new MyEdge(typename MyEdge::UnaryDaughters(this,rule));
  }


  assert(rule);
  assert(rule->get_probability().size() > 0);

  assert(rule);
  assert(rule->get_probability().size() > 0);
  (*e)->get_annotations().inside_probabilities_unary_temp.array[0] += L_inside * rule->get_probability()[0][0];
}

template<class MyEdge>
void PCKYAllCell<MyEdge>::reset_probabilities()
{
  //  std::cout << max_size << std::endl;
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      edges[i]->reset_probabilities(0.0);
    }
}


template <class MyEdge>
void PCKYAllCell<MyEdge>::adjust_inside_probability()
{
  apply_on_edges(&Edge::adjust_inside_probability);
}




template<class MyEdge>
void PCKYAllCell<MyEdge>::compute_inside_probabilities()
{
//   apply_on_edges( & Edge::clean_invalidated_binaries);
  
  apply_on_edges(std::function<void(Edge&)>([](Edge& edge){if (edge.get_lex()) edge.get_annotations().reset_probabilities();}) ,
                      & Edge::LexicalDaughters::update_inside_annotations  ,
                      & Edge:: BinaryDaughters::update_inside_annotations  ,
                      & Edge::                  prepare_inside_probability );
  
  apply_on_edges(& Edge::  UnaryDaughters::update_inside_annotations);
  apply_on_edges(& Edge::                  adjust_inside_probability);
}


template<class MyEdge>
void PCKYAllCell<MyEdge>::compute_outside_probabilities()
{
  apply_on_edges(& Edge::             prepare_outside_probability);
  apply_on_edges(& Edge::UnaryDaughters::update_outside_annotations);
  apply_on_edges(& Edge::              adjust_outside_probability);
  apply_on_edges(& Edge::BinaryDaughters::update_outside_annotations);
}


///////////



template<class MyEdge>
void PCKYAllCell<MyEdge>::compute_best_viterbi_derivation(const AnnotatedLabelsInfo& symbol_map)
{
  //iterate through all the packed edges in this cell, processing the lexical daughters first
  for(unsigned i = 0; i < max_size; ++i) {
    if(exists_edge(i))  {
      edges[i]->create_viterbi(symbol_map.get_number_of_annotations(i));
      edges[i]->replace_rule_probabilities(0);
      edges[i]->compute_best_lexical();
      edges[i]->compute_best_binary();
    }
  }
  //now process the unary daughters
  for(unsigned i = 0; i < max_size; ++i) {
    if(exists_edge(i)) {
      edges[i]->compute_best_unary();
    }
  }
}

// these 2 predicates returns true if the daughter(s) can be removed
// ie, if it's pointing to invalid edges
template <typename Cell>
struct pred_beam_clean_bin : public std::unary_function<typename Cell::Edge::BinaryDaughters, bool>
{
  bool operator()(const typename Cell::Edge::BinaryDaughters& packededgedaughter) const
  {

    Cell * cell0 = packededgedaughter.left_daughter();
    if(cell0->is_closed()) return true;
    typename Cell::Edge * lefty = cell0->get_edge_ptr(packededgedaughter.get_rule()->get_rhs0());
    if (lefty == NULL) return true;


    Cell * cell1 = packededgedaughter.right_daughter();
    if(cell1->is_closed()) return true;
    typename Cell::Edge * righty = cell1->get_edge_ptr(packededgedaughter.get_rule()->get_rhs1());
    return righty == NULL;
  }
};


template <typename Cell>
struct pred_beam_clean_un : public std::unary_function<typename Cell::Edge::UnaryDaughters, bool>
{
  bool operator()(const typename Cell::Edge::UnaryDaughters& packededgedaughter) const
  {
    Cell * cell = packededgedaughter.left_daughter();
    // cell should be equal to the current cell so this test is useless
    //    if(cell->is_closed()) return true;
    return cell->get_edge_ptr(packededgedaughter.get_rule()->get_rhs0()) == NULL;
  }
};

template<class MyEdge>
void PCKYAllCell<MyEdge>::clean()
{

  bool changed;
  do {
    changed =  false;

    // go through all the lists of unary daughters and remove the ones pointing on removed edges
    for(unsigned i = 0; i < max_size; ++i)
      if(edges[i]) {
	std::vector<typename MyEdge::UnaryDaughters >& udaughters = edges[i]->get_unary_daughters();
	udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(),
                                        pred_beam_clean_un<PCKYAllCell<MyEdge> >()),
			 udaughters.end());

	if(edges[i]->get_binary_daughters().empty() &&
           edges[i]->get_lexical_daughters().empty() &&
           edges[i]->get_unary_daughters().empty()) {
	  //	std::cout << "I shall be removed!" << std::endl;
	  delete edges[i];
	  edges[i]=NULL;
	  changed =  true;
	}
      }
  }
  while(changed);

  // final memory reclaim
  // TODO: benchmark this carefully
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      std::vector<typename MyEdge::UnaryDaughters >& udaughters = edges[i]->get_unary_daughters();
      if(udaughters.capacity() != udaughters.size()) {
        std::vector<typename MyEdge::UnaryDaughters > tmp;
        tmp.swap(udaughters);
        udaughters.insert(udaughters.begin(), tmp.begin(), tmp.end());
      }
    }

  // //should be a proper method
  // //if all edge pointers are NULL, close the cell
  //an overcomplicated way to do the same as below
  if(std::find_if(edges,edges+max_size,std::bind2nd(std::not_equal_to<MyEdge*>(), (MyEdge*)NULL)) == edges+max_size)
    {
      closed = true;
      delete[] edges;
      edges = NULL;
    }

  // bool all_null = true;
  // for(unsigned i = 0; i < max_size; ++i) {
  //   if (edges[i]) {
  //     all_null = false;
  //     break;
  //   }
  // }
  // if(all_null) {
  //   //  std::cout << "ALL NULL" << std::endl;
  //   closed = true;
  //   delete edges;
  //   edges = NULL;
  // }
}




//relative prior beam
template <class MyEdge>
void PCKYAllCell<MyEdge>::beam(const std::vector<double>& priors, double threshold)
{
  double max = 0.0;
  double beam = threshold;

  std::vector<double> sums = priors;

  //computing unannotated inside probabilities
  //looking for the probablity of the most probable symbol
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      sums[i] *= std::accumulate(edges[i]->get_annotations().inside_probabilities.array.begin(),
                                 edges[i]->get_annotations().inside_probabilities.array.end(),
                                 0.0);
      max = std::max(max, sums[i]);
    }

  //setting threshold
  beam *= max;

  //looking for edges below threshold
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      if(sums[i] < beam) {
        delete edges[i];
        edges[i]=NULL;
      }
    }

  //  clean the cell

  clean();
}



// Relative Inside/Outside beam
template<class MyEdge>
void PCKYAllCell<MyEdge>::beam(double threshold)
{
  double max = 0.0;
  double beam = threshold;

  std::vector<double> sums(max_size,0.0);

  //computing unannotated inside probabilities
  //looking for the probability of the most probable symbol
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      double ins = std::accumulate(edges[i]->get_annotations().inside_probabilities.array.begin(),
                                   edges[i]->get_annotations().inside_probabilities.array.end(),
                                   0.0);
      double outs = std::accumulate(edges[i]->get_annotations().outside_probabilities.array.begin(),
                                    edges[i]->get_annotations().outside_probabilities.array.end(),
                                    0.0);

      sums[i] = ins * outs;
      if(max < sums[i]) {max = sums[i];}
    }

  //setting threshold
  beam *= max;

  //looking for edges below threshold
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      if(sums[i] < beam) {
        delete edges[i];
        edges[i]=NULL;
      }
    }

  //  clean the cell
  clean();
}


template<class MyEdge>
void PCKYAllCell<MyEdge>::clean_binary_daughters()
{
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      MyEdge * edge = edges[i];

      // go through all the lists of binary daughters and remove the ones pointing on removed edges
      std::vector<typename MyEdge::BinaryDaughters >& bdaughters = edge->get_binary_daughters();

      bdaughters.erase(std::remove_if(bdaughters.begin(), bdaughters.end(), pred_beam_clean_bin<PCKYAllCell<MyEdge> >()), bdaughters.end());


      // Reclaim memory !
      if(bdaughters.capacity() != bdaughters.size()) {
        std::vector<typename MyEdge::BinaryDaughters > tmp;
        tmp.swap(bdaughters);
        bdaughters.insert(bdaughters.begin(), tmp.begin(), tmp.end());
      }
    }
}


// Absolute Inside/Outside beam
template<class MyEdge>
void PCKYAllCell<MyEdge>::beam(double log_threshold, double log_sent_prob)
{
  double beam = log_threshold  + log_sent_prob;

  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      bool all_invalid = true;
      AnnotationInfo& ai = edges[i]->get_annotations();

      // calculate posterior for each annotation
      for(unsigned annot = 0 ; annot < ai.inside_probabilities.array.size(); ++annot) {
        if(ai.inside_probabilities.array[annot] != LorgConstants::NullProba
          //|| ai.outside_probabilities.array[annot] != LorgConstants::NullProba
        ) {

          double prob = std::log(ai.inside_probabilities.array[annot]) + std::log(ai.outside_probabilities.array[annot]);
          //          double prob = std::log(ai.inside_probabilities.array[annot] * ai.outside_probabilities.array[annot]);

          if (prob > beam)
            all_invalid = false;
          else {
            ai.inside_probabilities.array[annot] = ai.outside_probabilities.array[annot] = LorgConstants::NullProba;
          }
        }
      }

      //remove edge if all annotations are NullProba
      if(all_invalid) {
        delete edges[i];
        edges[i]=NULL;
      }
    }
  // you must call clean after this method
}


// returns true if the branching can be removed
// in the sense of Huang, 2008
template <typename Cell>
struct pred_beam_huang
{
  double log_threshold;
  double log_outside_up;

  pred_beam_huang(double th, double se, double ou) : log_threshold(th + se), log_outside_up(ou) {}


  // assume that clean has already been called
  // and so lefty and righty are never NULL
  bool operator()(const typename Cell::Edge::BinaryDaughters& packededgedaughter) const
  {


    Cell * cell0 = packededgedaughter.left_daughter();
    assert(cell0 != NULL);

    typename Cell::Edge * lefty = cell0->get_edge_ptr(packededgedaughter.get_rule()->get_rhs0());
    assert(lefty != NULL);
    const AnnotationInfo& ailefty = lefty->get_annotations();

    double total_in = 0;
    double sum = 0;
    for(unsigned annot = 0 ; annot < ailefty.inside_probabilities.array.size(); ++annot) {
      if(ailefty.inside_probabilities.array[annot] != LorgConstants::NullProba) {
        sum += ailefty.inside_probabilities.array[annot];
      }
    }
    total_in += std::log(sum);


    Cell * cell1 = packededgedaughter.right_daughter();
    assert(cell1 != NULL);
    typename Cell::Edge * righty = cell1->get_edge_ptr(packededgedaughter.get_rule()->get_rhs1());
    assert(righty != NULL);
    const AnnotationInfo& airighty = righty->get_annotations();

    sum = 0;
    for(unsigned annot = 0 ; annot < airighty.inside_probabilities.array.size(); ++annot) {
      if(airighty.inside_probabilities.array[annot] != LorgConstants::NullProba) {
        sum += airighty.inside_probabilities.array[annot];
      }
    }

    total_in += std::log(sum);

    bool remove = log_outside_up + total_in  < log_threshold;

    return remove;
  }

  bool operator()(const typename Cell::Edge::UnaryDaughters& packededgedaughter) const
  {
    Cell * cell = packededgedaughter.left_daughter();
    assert(cell != NULL);
    typename Cell::Edge * lefty = cell->get_edge_ptr(packededgedaughter.get_rule()->get_rhs0());
    assert(lefty != NULL);

    const AnnotationInfo& ailefty = lefty->get_annotations();

    double total_in = 0;

    for(unsigned annot = 0 ; annot < ailefty.inside_probabilities.array.size(); ++annot) {
      if(ailefty.inside_probabilities.array[annot] != LorgConstants::NullProba) {
        total_in += ailefty.inside_probabilities.array[annot];
      }
    }

    total_in = std::log(total_in);

    bool remove = log_outside_up + total_in  < log_threshold;

    return remove;

  }
};

template<class MyEdge>
void PCKYAllCell<MyEdge>::beam_huang(double log_threshold, double log_sent_prob)
{
  for(unsigned i = 0; i < max_size; ++i) {
    // std::cout << edges << std::endl;
    // std::cout << i << std::endl;
    if(edges[i]) {
      MyEdge *  edge = edges[i];
      AnnotationInfo& ai = edge->get_annotations();

      double total_out = 0;
      for(unsigned annot = 0 ; annot < ai.outside_probabilities.array.size(); ++annot) {
        if(ai.outside_probabilities.array[annot] != LorgConstants::NullProba) {
          total_out += ai.outside_probabilities.array[annot];
        }
      }

      total_out = std::log(total_out);

      pred_beam_huang<PCKYAllCell<MyEdge> > huang(log_threshold, log_sent_prob, total_out);


      std::vector<typename MyEdge::BinaryDaughters >& bdaughters = edge->get_binary_daughters();
      bdaughters.erase(std::remove_if(bdaughters.begin(), bdaughters.end(), huang),
                       bdaughters.end());

      std::vector<typename MyEdge::UnaryDaughters >& udaughters = edge->get_unary_daughters();

      udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(), huang),
                       udaughters.end());
    }
  }
}

template<class MyEdge>
void PCKYAllCell<MyEdge>::change_rules_resize(const AnnotatedLabelsInfo& next_annotations,
                                              const std::vector<std::vector<std::vector<unsigned> > >& annot_descendants_current)
{
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {

      AnnotationInfo a(next_annotations.get_number_of_annotations(i), 0.0);

      //process invalid annotations
      for(unsigned annot = 0; annot < edges[i]->get_annotations().inside_probabilities.array.size(); ++annot) {
        if(!edges[i]->valid_prob_at(annot)) {

          const std::vector<unsigned>& next_invalids = annot_descendants_current[i][annot];
          for(std::vector<unsigned>::const_iterator new_annot(next_invalids.begin()); new_annot != next_invalids.end(); ++new_annot) {
            a.inside_probabilities.array[*new_annot] = LorgConstants::NullProba;
            a.outside_probabilities.array[*new_annot] = LorgConstants::NullProba;
          }


        }
      }

      //replace annot
      std::swap(a,edges[i]->get_annotations());

      //replace rule
      edges[i]->replace_rule_probabilities(0);

    }
}



template<class MyEdge>
void PCKYAllCell<MyEdge>::change_rules_resize(unsigned new_size, unsigned finer_idx)
{
  for(unsigned i = 0; i < max_size; ++i)
    {
      if(edges[i]) {
        //resize
        edges[i]->get_annotations().reset_probabilities(0.0);
        edges[i]->get_annotations().resize(new_size);


        //replace rule
        edges[i]->replace_rule_probabilities(finer_idx);
      }
    }
}



//simple stuff
template<class MyEdge>
std::ostream& operator<<(std::ostream& out, const PCKYAllCell<MyEdge>& cell)
{
  int nb_entries = 0;
  for(unsigned i = 0; i < cell.max_size ; ++i)
    if(cell.edges[i]) {
      ++nb_entries;
      out << i << ":" << *cell.edges[i] << std::endl;
    }
  return out << "filled entries: " << nb_entries;
}

template<class MyEdge>
void PCKYAllCell<MyEdge>::clear()
{
  closed = false;

  if(!edges)
    {
      edges =  new MyEdge * [max_size];
      memset(edges, 0, max_size * sizeof(MyEdge*));
    }
  else
    for(unsigned i = 0; i < max_size; ++i)
      {
        if(edges[i]) {
          delete edges[i];
          edges[i] = 0;
        }
      }
}



#endif //PCKYALLCELL_H
