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
void PCKYAllCell<Types>::clear()
{
  apply_on_edges( & Edge::close );
}





/* *****************************************************************************************************************
 *                                             ITERATEURS                                                          *
 * *****************************************************************************************************************/


template<class Edge>
Edge * begin(Edge * e) { return e; }

template<class Edge>
Edge * end(Edge * e) { return (e+Edge::Cell::get_max_size()); }

template<class Edge, class PEdge>
struct _pedge_iterator {
  Edge * it ;
  PEdge & operator*() { return *(PEdge *)it; }
  _pedge_iterator<Edge,PEdge> & operator++() { ++it; return *this; }
  _pedge_iterator<Edge,PEdge>(PEdge & init) { it = (Edge *) (&init); }
  bool operator!=(const _pedge_iterator< Edge, PEdge > &o) const { return it!=o.it; }
};

template<class Edge>
struct Uedges {
  typename Edge::UEdge * beg;
  typename Edge::UEdge * en;
  _pedge_iterator<Edge,typename Edge::UEdge> & begin() { return _pedge_iterator<Edge,typename Edge::UEdge>(*beg); }
  _pedge_iterator<Edge,typename Edge::UEdge> &   end() { return _pedge_iterator<Edge,typename Edge::UEdge>(*en); }
  Uedges(Edge * e) : beg(&e->uedge()), en(&(e+Edge::Cell::get_max_size())->uedge()) {}
};

template<class Edge>
struct LBedges {
  typename Edge::LBEdge * beg;
  typename Edge::LBEdge * en;
  _pedge_iterator<Edge,typename Edge::UEdge> & begin() { return _pedge_iterator<Edge,typename Edge::LBEdge>(*beg); }
  _pedge_iterator<Edge,typename Edge::UEdge> &   end() { return _pedge_iterator<Edge,typename Edge::LBEdge>(*en); }
  LBedges(Edge * e) : beg(&e->lbedge()), en(&(e+Edge::Cell::get_max_size())->lbedge()) {}
};








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
                                           double R_inside)
{
  Edge & e = edges[rule->get_lhs()];
//   std::clog << "PCKYAllCell<Types>::process_candidate add_daughters" << &e << ", " <<&left <<", "<<&right <<")" <<std::endl;
  e.add_daughters(left,right,rule);
//   std::clog << "PCKYAllCell<Types>::process_candidate inside_probabilities" << &e << ", " <<&left <<", "<<&right <<")" <<std::endl;
  e.get_annotations().inside_probabilities.array[0] += R_inside * rule->get_probability()[0][0][0];
}

template<class Types>
void PCKYAllCell<Types>::process_candidate(const UnaryRule* rule, double R_inside)
{
  assert(rule);
  assert(rule->get_probability().size() > 0);
  static int i = 0;
  ++i;

  Edge & e = edges[rule->get_lhs()];
  e.add_daughters(edges[rule->get_rhs0()],rule);

  //   std::cout << "PCKYAllCell<Types>::process_candidate. array at " << & e.get_annotations().inside_probabilities_unary_temp.array[0] << std::endl; std::cout.flush();
  e.uedge().get_annotations().inside_probabilities.array[0] += R_inside * rule->get_probability()[0][0];
}

template<class Types>
inline
void PCKYAllCell<Types>::process_candidates(const Word & word)
{
//   assert(edges.size() == max_size);
  for(const auto & rule : word.get_rules())
  {
    const typename Types::LRule* r = static_cast<const typename Types::LRule*>(rule);
    Edge & e = edges[rule->get_lhs()];
    e.add_daughters(r,&word);
    e.get_annotations().inside_probabilities.array[0] += r->get_probability()[0];
  }
}


template<class Types>
void PCKYAllCell<Types>::reset_probabilities()
{
  apply_on_edges(toFunc(&Edge::reset_probabilities));
}





template<class Types>
void PCKYAllCell<Types>::compute_merged_inside_probabilities()
{
#warning usefull ?
// warning comment : one of the two following lines is useless. Maybe both.
// this is linked with the c_r_s function in PCKYAllCell<Types>::change_rules_resize
// which is called by ParserCKYAll_Impl<Types>::change_rules_resize
// which is called at the end of ParserCKYAll_Impl<Types>::beam_c2f
// The c_r_s function extends the annotation probabilities vectors, and nullify them,
// excepted for LorgConstants::NullProba, that are inherited from the previous pass in beam_c2f

//   apply_on_edges(toFunc(&Edge::reset_probabilities)); // useful 1 ?
  apply_on_lbedges(std::function<void(LBEdge&)>([](LBEdge& edge){if (edge.get_lex()) edge.get_annotations().reset_probabilities();})); // usefull 2 ?

  apply_on_lbedges(& Edge::update_merged_inside_annotations_from_lex  ,
                   & Edge::update_merged_inside_annotations_from_bin);

  apply_on_uedges(& Edge::update_unary_inside_annotations,
                  & Edge::add_unary_insides_to_merged  );
}

template<class Types>
void PCKYAllCell<Types>::compute_merged_outside_probabilities()
{
  apply_on_uedges( & Edge::copy_merged_outsides_to_unary,
                   & Edge::update_unary_outside_annotations);
  apply_on_lbedges(& Edge::update_binary_outside_annotations_from_merged);
}



template<class Types>
void PCKYAllCell<Types>::compute_inside_probabilities()
{
  apply_on_edges(toFunc(&Edge::reset_probabilities)); // useful 1 ?
  apply_on_lbedges(& Edge::update_lexical_inside_annotations,
                   & Edge:: update_binary_inside_annotations);
  apply_on_uedges (& Edge::  update_unary_inside_annotations);
}

template<class Types>
void PCKYAllCell<Types>::compute_outside_probabilities()
{
  apply_on_uedges( & Edge:: update_unary_outside_annotations);
  apply_on_lbedges(& Edge::update_binary_outside_annotations);
}

///////////



template<class Types>
void PCKYAllCell<Types>::clean()
{

  // go through all the lists of unary daughters and remove the ones pointing on removed edges
  for(auto & edge : edges)
    if(edge.uedge().is_opened()) {
      auto & udaughters = edge.uedge().get_unary_daughters();
      udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(),
                                      toFunc(& UnaryDaughter::points_towards_invalid_edges)),
                       udaughters.end());
      
      if (edge.uedge().no_daughters())
      {
        BLOCKTIMING("PCKYAllCell<Types>::clean - edge.close_u()");
        edge.close_u();
        if (edge.is_closed()) {
          BLOCKTIMING("PCKYAllCell<Types>::clean - edge.close_u() caused edge.is_closed()");
        }
      }
    }

  // final memory reclaim
  // TODO: benchmark this carefully
  bool all_null = true;
  for(auto & edge : edges)
    if(edge.uedge().is_opened()) {
      edge.uedge().get_unary_daughters().shrink_to_fit();
      all_null = false ;
    }
    else if(edge.is_opened())
      all_null = false ;

  // //if all edges are closed, close the cell
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
    if(edges[i].is_opened()) {
      sums[i] *= std::accumulate(edges[i].get_annotations().inside_probabilities.array.begin(),
                                 edges[i].get_annotations().inside_probabilities.array.end(),
                                 0.0);
      max = std::max(max, sums[i]);
    }

  //setting threshold
  beam *= max;
//   ubeam *= umax;

  //looking for edges below threshold
  for(unsigned i = 0; i < max_size; ++i) {
    if(edges[i].is_opened()) {
      if(sums[i] < beam) {
        BLOCKTIMING("PCKYAllCell<Types>::beam(priors) - edges[i].close");
        edges[i].close();
      }
    }
  }
  //  clean the cell

  clean();
}



template<class EdgeType>
static void absolute_beam(EdgeType & edge, double beam)
{
  BLOCKTIMING("absolute_beam");
   bool all_invalid = true;
    AnnotationInfo& ai = edge.get_annotations();

    // calculate posterior for each annotation
    for(unsigned annot = 0 ; annot < ai.inside_probabilities.array.size(); ++annot) {
      if(ai.inside_probabilities.array[annot] != LorgConstants::NullProba
        //|| ai.outside_probabilities.array[annot] != LorgConstants::NullProba
      ) {

        double prob = std::log(ai.inside_probabilities.array[annot]) + std::log(ai.outside_probabilities.array[annot]);
        //          double prob = std::log(ai.inside_probabilities.array[annot] * ai.outside_probabilities.array[annot]);

        if (prob > beam) {
          BLOCKTIMING("absolute_beam --> remain");
          all_invalid = false;
        }
        else {
          BLOCKTIMING("absolute_beam --> NullProba");
          ai.inside_probabilities.array[annot] = ai.outside_probabilities.array[annot] = LorgConstants::NullProba;
        }
      }
    }

    //remove edge if all annotations are NullProba
    if(all_invalid) {
      BLOCKTIMING("absolute_beam --> edge.close");
     edge.close();
    }
}

// Absolute Inside/Outside beam
template<class Types>
void PCKYAllCell<Types>::beam(double log_threshold, double log_sent_prob)
{
  double beam = log_threshold  + log_sent_prob;
  apply_on_edges(function<void(Edge&)>(std::bind(absolute_beam<Edge>, std::placeholders::_1, beam)));
}


// returns true if the branching can be removed
// in the sense of Huang, 2008
template <typename Edge>
struct pred_beam_huang_unary
{
  double log_threshold;
  double log_outside_up;

  pred_beam_huang_unary(double th, double se, double ou) : log_threshold(th + se), log_outside_up(ou) {}

  bool operator()(const typename Edge::UnaryDaughter& packededgedaughter) const
  {
    Edge & lefty = packededgedaughter.daughter();
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

template <typename Edge>
struct pred_beam_huang_binary
{
  double log_threshold;
  double log_outside_up;

  pred_beam_huang_binary(double th, double se, double ou) : log_threshold(th + se), log_outside_up(ou) {}


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
};


template<class Types>
void PCKYAllCell<Types>::beam_huang(double log_threshold, double log_sent_prob)
{
  for(Edge & edge: edges) {
    // std::cout << edges << std::endl;
    // std::cout << i << std::endl;
    if(not edge.lbedge().is_closed()) {
      AnnotationInfo& ai = edge.lbedge().get_annotations();

      double total_out = 0;
      for(unsigned annot = 0 ; annot < ai.outside_probabilities.array.size(); ++annot) {
        if(ai.outside_probabilities.array[annot] != LorgConstants::NullProba) {
          total_out += ai.outside_probabilities.array[annot];
        }
      }

      total_out = std::log(total_out);

      pred_beam_huang_binary<typename Types::AEdge> huang(log_threshold, log_sent_prob, total_out);

      std::vector<typename Types::BinaryDaughter >& bdaughters = edge.lbedge().get_binary_daughters();
      bdaughters.erase(std::remove_if(bdaughters.begin(), bdaughters.end(), huang),
                       bdaughters.end());
    }
    if(not edge.uedge().is_closed()) {
      AnnotationInfo& ai = edge.uedge().get_annotations();

      double total_out = 0;
      for(unsigned annot = 0 ; annot < ai.outside_probabilities.array.size(); ++annot) {
        if(ai.outside_probabilities.array[annot] != LorgConstants::NullProba) {
          total_out += ai.outside_probabilities.array[annot];
        }
      }

      total_out = std::log(total_out);

      pred_beam_huang_unary<typename Types::AEdge> huang(log_threshold, log_sent_prob, total_out);

      std::vector<typename Types::UnaryDaughter >& udaughters = edge.uedge().get_unary_daughters();

      udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(), huang),
                       udaughters.end());
    }
  }
}

template<class Types>
void PCKYAllCell<Types>::change_rules_resize(const AnnotatedLabelsInfo& next_annotations,
                                              const std::vector<std::vector<std::vector<unsigned> > >& annot_descendants_current)
{
  auto c_r_s = [&](AEdge & edge, int i)->bool {
    if(not edge.is_closed()) {
      
      AnnotationInfo a(next_annotations.get_number_of_annotations(i), 0.0);

      //process invalid annotations
      for(unsigned annot = 0; annot < edge.get_annotations().inside_probabilities.array.size(); ++annot) {
        if(!edge.valid_prob_at(annot)) {

          for(auto new_annot : annot_descendants_current[i][annot]) {
            a.inside_probabilities.array[new_annot] = 
            a.outside_probabilities.array[new_annot] = LorgConstants::NullProba;
          }
        }
      }

      //replace annot
      edge.get_annotations() = std::move(a);
      return true;
    }
    return false;
  };
  
  for(size_t i=0; i<max_size; ++i) {
    Edge & edge = edges[i];
    if (c_r_s(edge,          i)) {
      if (c_r_s(edge. uedge(), i)) edge. uedge().replace_rule_probabilities(0);
      if (c_r_s(edge.lbedge(), i)) edge.lbedge().replace_rule_probabilities(0);
    }
  }
}



template<class Types>
void PCKYAllCell<Types>::change_rules_resize(unsigned new_size, unsigned finer_idx)
{
  // resize annotations
  apply_on_edges  ( function<void(Edge&)>([new_size](Edge&e){ e.resize_annotations(new_size); }));
  //replace rules
  apply_on_lbedges( function<void(LBEdge&)>([finer_idx](LBEdge&e){
    e.replace_rule_probabilities(finer_idx);
  })
  );
  apply_on_uedges( function<void(UEdge&)>([finer_idx](UEdge&e){
    e.replace_rule_probabilities(finer_idx);
  })
  );
}



//simple stuff
template<class Types>
void PCKYAllCell<Types>::dump(std::ostream& out) const
{
  out << "(cell: span=" << get_end() - get_begin() + 1 << ", beg=" << get_begin() << " :"<< std::endl;
  int nb_entries = 0;
  for(unsigned i = 0; i < max_size ; ++i)
    if(not edges[i].is_closed()) {
      ++nb_entries;
      out << " " << i << ":" << edges[i] << std::endl;
    }
  out << "filled entries: " << nb_entries << ")";
}





#endif //PCKYALLCELL_HPP
