// -*- mode: c++ -*-
#ifndef PCKYALLCELL_HPP
#define PCKYALLCELL_HPP

#include "PCKYAllCell.h"
#include "edges/AnnotationInfo.h"

#include <cassert>
#include <cstring>

#include <numeric>
#include <algorithm>
#include <functional>


using std::function;


template<class Types>
unsigned PCKYAllCell<Types>::max_size = 0;


template<class Types>
PCKYAllCell<Types>::~PCKYAllCell()
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

template<class Types>
void PCKYAllCell<Types>::process_candidate(Cell * left,
                                            Cell * right,
                                            const BinaryRule* rule,
                                            double LR_inside)
{
  Edge ** e = &edges[rule->get_lhs()];

  if(*e)
    (*e)->add_daughters(left,right,rule);
  else {
    *e = new Edge(BinaryDaughter(left,right,rule));
  }


  (*e)->get_annotations().inside_probabilities.array[0] += LR_inside * rule->get_probability()[0][0][0];
}

template<class Types>
void PCKYAllCell<Types>::process_candidate(const UnaryRule* rule, double L_inside)
{
  assert(rule);
  assert(rule->get_probability().size() > 0);


  typename Types::Edge ** e = &edges[rule->get_lhs()];

  if(*e)  {
    (*e)->add_daughters(this,rule);
  }
  else {
    //std::cout <<" adding a new edge " << *rule << std::endl;
    *e = new Edge(UnaryDaughter(this,rule));
  }


  assert(rule);
  assert(rule->get_probability().size() > 0);

  assert(rule);
  assert(rule->get_probability().size() > 0);
  (*e)->get_annotations().inside_probabilities_unary_temp.array[0] += L_inside * rule->get_probability()[0][0];
}

template<class Types>
void PCKYAllCell<Types>::reset_probabilities()
{
  apply_on_edges(function<void(Edge&)>([](Edge&e){e.get_annotations().reset_probabilities(0.0);}));
}


template <class Types>
void PCKYAllCell<Types>::adjust_inside_probability()
{
  apply_on_edges(&Edge::adjust_inside_probability);
}




template<class Types>
void PCKYAllCell<Types>::compute_inside_probabilities()
{
  //   apply_on_edges( & Edge::clean_invalidated_binaries);

  apply_on_edges(std::function<void(Edge&)>([](Edge& edge){if (edge.get_lex()) edge.get_annotations().reset_probabilities();}) ,
                      & LexicalDaughter::update_inside_annotations  ,
                      &  BinaryDaughter::update_inside_annotations  ,
                      &            Edge::prepare_inside_probability );

  apply_on_edges(& UnaryDaughter::update_inside_annotations);
  apply_on_edges(& Edge::         adjust_inside_probability);
}

template<class Types>
void PCKYAllCell<Types>::compute_outside_probabilities()
{
  apply_on_edges(& Edge::             prepare_outside_probability);
  apply_on_edges(&     UnaryDaughter::update_outside_annotations);
  apply_on_edges(& Edge::              adjust_outside_probability);
  apply_on_edges(&    BinaryDaughter::update_outside_annotations);
}


///////////



template<class Types>
void PCKYAllCell<Types>::clean()
{

  bool changed;
  do {
    changed =  false;

    // go through all the lists of unary daughters and remove the ones pointing on removed edges
    for(unsigned i = 0; i < max_size; ++i)
      if(edges[i]) {
        auto & udaughters = edges[i]->get_unary_daughters();
        udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(),
                                        toFunc(& UnaryDaughter::points_towards_invalid_cells)),
                         udaughters.end());

        if (edges[i]->get_binary_daughters().empty()
            && edges[i]->get_lexical_daughters().empty()
            && edges[i]->get_unary_daughters().empty())
        {
          //std::cout << "I shall be removed!" << std::endl;
          delete edges[i];
          edges[i]=NULL;
          changed =  true;
        }
      }
  } while(changed);

  // final memory reclaim
  // TODO: benchmark this carefully
  bool all_null = true;
  for(unsigned i = 0; i < max_size; ++i)
    if(edges[i]) {
      std::vector<UnaryDaughter >& udaughters = edges[i]->get_unary_daughters();
      if(udaughters.capacity() != udaughters.size()) {
        std::vector<UnaryDaughter> tmp;
        tmp.swap(udaughters);
        udaughters.insert(udaughters.begin(), tmp.begin(), tmp.end());
      }
      all_null = false ;
    }

  // //if all edge pointers are NULL, close the cell
  if(all_null)
    {
      closed = true;
      delete[] edges;
      edges = NULL;
    }
}




//relative prior beam
template <class Types>
void PCKYAllCell<Types>::beam(const std::vector<double>& priors, double threshold)
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
template<class Types>
void PCKYAllCell<Types>::beam(double threshold)
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


// Absolute Inside/Outside beam
template<class Types>
void PCKYAllCell<Types>::beam(double log_threshold, double log_sent_prob)
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
  bool operator()(const typename Cell::Edge::BinaryDaughter& packededgedaughter) const
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

  bool operator()(const typename Cell::Edge::UnaryDaughter& packededgedaughter) const
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

template<class Types>
void PCKYAllCell<Types>::beam_huang(double log_threshold, double log_sent_prob)
{
  for(unsigned i = 0; i < max_size; ++i) {
    // std::cout << edges << std::endl;
    // std::cout << i << std::endl;
    if(edges[i]) {
      typename Types::Edge *  edge = edges[i];
      AnnotationInfo& ai = edge->get_annotations();

      double total_out = 0;
      for(unsigned annot = 0 ; annot < ai.outside_probabilities.array.size(); ++annot) {
        if(ai.outside_probabilities.array[annot] != LorgConstants::NullProba) {
          total_out += ai.outside_probabilities.array[annot];
        }
      }

      total_out = std::log(total_out);

      pred_beam_huang<PCKYAllCell<Types> > huang(log_threshold, log_sent_prob, total_out);


      std::vector<typename Types::BinaryDaughter >& bdaughters = edge->get_binary_daughters();
      bdaughters.erase(std::remove_if(bdaughters.begin(), bdaughters.end(), huang),
                       bdaughters.end());

      std::vector<typename Types::UnaryDaughter >& udaughters = edge->get_unary_daughters();

      udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(), huang),
                       udaughters.end());
    }
  }
}

template<class Types>
void PCKYAllCell<Types>::change_rules_resize(const AnnotatedLabelsInfo& next_annotations,
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



template<class Types>
void PCKYAllCell<Types>::change_rules_resize(unsigned new_size, unsigned finer_idx)
{
  apply_on_edges(function<void(Edge&)>([new_size,finer_idx](Edge&e){
    //resize
    e.get_annotations().reset_probabilities(0.0);
    e.get_annotations().resize(new_size);
    //replace rule
    e.replace_rule_probabilities(finer_idx);
  })
  );
}



//simple stuff
template<class Types>
std::ostream& operator<<(std::ostream& out, const PCKYAllCell<Types>& cell)
{
  int nb_entries = 0;
  for(unsigned i = 0; i < cell.max_size ; ++i)
    if(cell.edges[i]) {
      ++nb_entries;
      out << i << ":" << *cell.edges[i] << std::endl;
    }
  return out << "filled entries: " << nb_entries;
}

template<class Types>
void PCKYAllCell<Types>::clear()
{
  closed = false;

  if(!edges)
  {
    edges =  new typename Types::Edge * [max_size];
    memset(edges, 0, max_size * sizeof(typename Types::Edge*));
  }
  else
    for(unsigned i = 0; i < max_size; ++i)
    {
      if(edges[i]) {
        delete edges[i];
        edges[i] = nullptr;
      }
    }
}



#endif //PCKYALLCELL_HPP
