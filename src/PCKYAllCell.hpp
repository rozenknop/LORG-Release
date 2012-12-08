// -*- mode: c++ -*-
#ifndef PCKYALLCELL_HPP
#define PCKYALLCELL_HPP

#include "./PCKYAllCell.h"
#include "edges/AnnotationInfo.h"
#include "edges/PackedEdge.hpp"

#include <cassert>
#include <cstring>

#include <numeric>
#include <algorithm>
#include <functional>


using std::function;


template<class Types>
unsigned PCKYAllCell<Types>::max_size = 0;

template<class Types>
PCKYAllCell<Types>::PCKYAllCell() {
}

template<class Types>
PCKYAllCell<Types>::PCKYAllCell( const PCKYAllCell<Types> & o )
{
  edges = (Edge*) new char[max_size*sizeof(Edge)];
  memcpy(edges, o.edges, max_size*sizeof(Edge));

//   std::copy(&(o.edges), &(o.edges)+1, &edges);
//   std::cout << "size of edges = " << edges.size() << " " << edges.capacity() <<  " " << max_size << std:: endl;
//   memcpy(e, o.edges.data(), max_size*sizeof(Edge));
  
//   std::cout << "copy constructor of " << this << " from " << &o << std::endl;
//   *this = o ;
//   memcpy(this, &o, sizeof(PCKYAllCell<Types>));

}

template<class Types>
void PCKYAllCell<Types>::clear()
{
//   static Edge protoEdge;
//   edges.assign(max_size, protoEdge);
  std::fill((char*)edges, (char*) (edges+max_size), 0);
}




template<class Types>
PCKYAllCell<Types>::~PCKYAllCell()
{
  apply_on_edges( &Edge::close );
}

template<class Types>
void PCKYAllCell<Types>::reserve_binary_daughters(const std::vector<int> & counts)
{
  for(int i=counts.size()-1; i>=0; --i) {
    if (counts[i]!=0) {
      edges[i].reserve_binary_daughters(counts[i]);
    }
  }
}

template<class Types>
void PCKYAllCell<Types>::process_candidate(Edge & left,
                                            Edge & right,
                                            const BinaryRule* rule,
                                            double LR_inside)
{
  Edge & e = edges[rule->get_lhs()];
  e.add_daughters(left,right,rule);
  e.get_annotations().inside_probabilities.array[0] += LR_inside * rule->get_probability()[0][0][0];
}

template<class Types>
void PCKYAllCell<Types>::process_candidate(const UnaryRule* rule, double L_inside)
{
  assert(rule);
  assert(rule->get_probability().size() > 0);
  static int i = 0;
  ++i;

  Edge & e = edges[rule->get_lhs()];
  e.add_daughters(edges[rule->get_rhs0()],rule);

  //   std::cout << "PCKYAllCell<Types>::process_candidate. array at " << & e.get_annotations().inside_probabilities_unary_temp.array[0] << std::endl; std::cout.flush();
  e.get_annotations().inside_probabilities_unary_temp.array[0] += L_inside * rule->get_probability()[0][0];
  
}

template<class Types>
inline
void PCKYAllCell<Types>::add_word(const Word & word)
{
  typedef typename Types::LexicalDaughter LDaughters;
//   assert(edges.size() == max_size);
  for(const auto & rule : word.get_rules())
  {
    const typename Types::LRule* r = static_cast<const typename Types::LRule*>(rule);
    int tag = rule->get_lhs();
    if (0==edges[tag].get_annotations().get_size())
      edges[tag].local_resize_annotations(1);
    edges[tag].add_daughters(r, &word);
    edges[tag].get_annotations().inside_probabilities.array[0] += r->get_probability()[0];
  }
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
//     apply_on_edges( & Edge::clean_invalidated_binaries);

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
    for(auto & edge : edges)
      if(not edge.is_closed()) {
        auto & udaughters = edge.get_unary_daughters();
        udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(),
                                        toFunc(& UnaryDaughter::points_towards_invalid_edges)),
                         udaughters.end());

        if (edge.no_daughters())
        {
          edge.close();
          changed =  true;
        }
      }
  } while(changed);

  // final memory reclaim
  // TODO: benchmark this carefully
  bool all_null = true;
  for(auto & edge : edges)
    if(not edge.is_closed()) {
      edge.get_unary_daughters().shrink_to_fit();
      all_null = false ;
    }

  // //if all edge pointers are NULL, close the cell
  if(all_null)
    {
      closed = true;
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
    if(not edges[i].is_closed()) {
      sums[i] *= std::accumulate(edges[i].get_annotations().inside_probabilities.array.begin(),
                                 edges[i].get_annotations().inside_probabilities.array.end(),
                                 0.0);
      max = std::max(max, sums[i]);
    }

  //setting threshold
  beam *= max;

  //looking for edges below threshold
  for(unsigned i = 0; i < max_size; ++i)
    if(not edges[i].is_closed()) {
      if(sums[i] < beam) {
        edges[i].close();
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

  for(Edge & edge: edges)
    if(not edge.is_closed()) {
      bool all_invalid = true;
      AnnotationInfo& ai = edge.get_annotations();

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
        edge.close();
      }
    }
  // you must call clean after this method
}


// returns true if the branching can be removed
// in the sense of Huang, 2008
template <typename Edge>
struct pred_beam_huang
{
  double log_threshold;
  double log_outside_up;

  pred_beam_huang(double th, double se, double ou) : log_threshold(th + se), log_outside_up(ou) {}


  // assume that clean has already been called
  // and so lefty and righty are never NULL
  bool operator()(const typename Edge::BinaryDaughter& packededgedaughter) const
  {


    const Edge & lefty = packededgedaughter.left_daughter();
    assert(not lefty.is_closed());
    const AnnotationInfo& ailefty = lefty.get_annotations();

    double total_in = 0;
    double sum = 0;
    for(unsigned annot = 0 ; annot < ailefty.inside_probabilities.array.size(); ++annot) {
      if(ailefty.inside_probabilities.array[annot] != LorgConstants::NullProba) {
        sum += ailefty.inside_probabilities.array[annot];
      }
    }
    total_in += std::log(sum);


    const Edge & righty = packededgedaughter.right_daughter();
    assert(not righty.is_closed());
    const AnnotationInfo& airighty = righty.get_annotations();

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

  bool operator()(const typename Edge::UnaryDaughter& packededgedaughter) const
  {
    Edge & lefty = packededgedaughter.left_daughter();
    assert(not lefty.is_closed());
    const AnnotationInfo& ailefty = lefty.get_annotations();

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
  for(Edge & edge: edges) {
    // std::cout << edges << std::endl;
    // std::cout << i << std::endl;
    if(not edge.is_closed()) {
      AnnotationInfo& ai = edge.get_annotations();

      double total_out = 0;
      for(unsigned annot = 0 ; annot < ai.outside_probabilities.array.size(); ++annot) {
        if(ai.outside_probabilities.array[annot] != LorgConstants::NullProba) {
          total_out += ai.outside_probabilities.array[annot];
        }
      }

      total_out = std::log(total_out);

      pred_beam_huang<typename Types::Edge> huang(log_threshold, log_sent_prob, total_out);


      std::vector<typename Types::BinaryDaughter >& bdaughters = edge.get_binary_daughters();
      bdaughters.erase(std::remove_if(bdaughters.begin(), bdaughters.end(), huang),
                       bdaughters.end());

      std::vector<typename Types::UnaryDaughter >& udaughters = edge.get_unary_daughters();

      udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(), huang),
                       udaughters.end());
    }
  }
}

template<class Types>
void PCKYAllCell<Types>::change_rules_resize(const AnnotatedLabelsInfo& next_annotations,
                                              const std::vector<std::vector<std::vector<unsigned> > >& annot_descendants_current)
{
  for(size_t i=0; i<max_size; ++i) {
    Edge & edge = edges[i];
    if(not edge.is_closed()) {
      
      AnnotationInfo a(next_annotations.get_number_of_annotations(i), 0.0);

      //process invalid annotations
      for(unsigned annot = 0; annot < edge.get_annotations().inside_probabilities.array.size(); ++annot) {
        if(!edge.valid_prob_at(annot)) {

          const std::vector<unsigned>& next_invalids = annot_descendants_current[i][annot];
          for(std::vector<unsigned>::const_iterator new_annot(next_invalids.begin()); new_annot != next_invalids.end(); ++new_annot) {
            a.inside_probabilities.array[*new_annot] = LorgConstants::NullProba;
            a.outside_probabilities.array[*new_annot] = LorgConstants::NullProba;
          }
        }
      }

      //replace annot
      std::swap(a,edge.get_annotations());

      //replace rule
      edge.replace_rule_probabilities(0);

    }
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
  out << "(cell: span=" << cell.get_end() - cell.get_begin() << ", beg=" << cell.get_begin() << " :"<< std::endl;
  int nb_entries = 0;
  for(unsigned i = 0; i < cell.max_size ; ++i)
    if(not cell.edges[i].is_closed()) {
      ++nb_entries;
      out << " " << i << ":" << cell.edges[i] << std::endl;
    }
  return out << "filled entries: " << nb_entries << ")";
}


#endif //PCKYALLCELL_HPP
