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
//   static Edge protoEdge;
//   edges.assign(max_size, protoEdge);
  std::fill((char*)edges, (char*) (edges+max_size), 0);
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
void PCKYAllCell<Types>::process_candidate(PEdge & left,
                                            PEdge & right,
                                            const BinaryRule* rule,
                                            double LR_inside)
{
  LBEdge & e = edges[rule->get_lhs()].lbedge();
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

  UEdge & e = edges[rule->get_lhs()].uedge();
  e.add_daughters(edges[rule->get_rhs0()].lbedge(),rule);

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
    if (0==edges[tag].lbedge().get_annotations().get_size())
      edges[tag].lbedge().local_resize_annotations(1);
    edges[tag].lbedge().add_daughters(r, &word);
    edges[tag].lbedge().get_annotations().inside_probabilities.array[0] += r->get_probability()[0];
  }
}


template<class Types>
void PCKYAllCell<Types>::reset_probabilities()
{
  apply_on_uedges(function<void(UEdge&)>([](UEdge&e){e.get_annotations().reset_probabilities(0.0);}));
  apply_on_lbedges(function<void(LBEdge&)>([](LBEdge&e){e.get_annotations().reset_probabilities(0.0);}));
}





template<class Types>
void PCKYAllCell<Types>::compute_inside_probabilities()
{
//     apply_on_edges( & Edge::clean_invalidated_binaries);

  apply_on_edges(std::function<void(LBEdge&)>([](LBEdge& edge){if (edge.get_lex()) edge.get_annotations().reset_probabilities();}) ,
                 & LexicalDaughter::update_inside_annotations  ,
                 &  BinaryDaughter::update_inside_annotations  );

  apply_on_edges(& UnaryDaughter::update_inside_annotations);
}

template<class Types>
void PCKYAllCell<Types>::compute_outside_probabilities()
{
  apply_on_edges(&     UnaryDaughter::update_outside_annotations);
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
      if(not edge.uedge().is_closed()) {
        auto & udaughters = edge.uedge().get_unary_daughters();
        udaughters.erase(std::remove_if(udaughters.begin(), udaughters.end(),
                                        toFunc(& UnaryDaughter::points_towards_invalid_edges)),
                         udaughters.end());

        if (edge.uedge().no_daughters())
        {
          edge.uedge().close();
          changed =  true;
        }
      }
  } while(changed);

  // final memory reclaim
  // TODO: benchmark this carefully
  bool all_null = true;
  for(auto & edge : edges)
    if(not edge.is_closed()) {
      edge.uedge().get_unary_daughters().shrink_to_fit();
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
  double lbmax, umax = lbmax = 0.0;
  double lbbeam, ubeam = lbbeam = threshold;

  std::vector<double> lbsums, usums = lbsums = priors;

  //computing unannotated inside probabilities
  //looking for the probablity of the most probable symbol
  for(unsigned i = 0; i < max_size; ++i)
    if(not edges[i].uedge().is_closed()) {
      usums[i] *= std::accumulate(edges[i].uedge().get_annotations().inside_probabilities.array.begin(),
                                  edges[i].uedge().get_annotations().inside_probabilities.array.end(),
                                  0.0);
      umax = std::max(umax, usums[i]);
      lbsums[i] *= std::accumulate(edges[i].lbedge().get_annotations().inside_probabilities.array.begin(),
                                   edges[i].lbedge().get_annotations().inside_probabilities.array.end(),
                                   0.0);
      lbmax = std::max(lbmax, lbsums[i]);
    }

  //setting threshold
  lbbeam *= lbmax;
  ubeam *= umax;

  //looking for edges below threshold
  for(unsigned i = 0; i < max_size; ++i) {
    if(not edges[i].uedge().is_closed()) {
      if(usums[i] < ubeam) {
        edges[i].uedge().close();
      }
    }
    if(not edges[i].lbedge().is_closed()) {
      if(lbsums[i] < lbbeam) {
        edges[i].lbedge().close();
      }
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

template<class EdgeType>
static void absolute_beam(EdgeType & edge, double beam)
{
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

// Absolute Inside/Outside beam
template<class Types>
void PCKYAllCell<Types>::beam(double log_threshold, double log_sent_prob)
{
  double beam = log_threshold  + log_sent_prob;
  apply_on_uedges(function<void(UEdge&)>(std::bind(absolute_beam<UEdge>, std::placeholders::_1, beam)));
  apply_on_lbedges(function<void(LBEdge&)>(std::bind(absolute_beam<LBEdge>, std::placeholders::_1, beam)));
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

      pred_beam_huang_binary<typename Types::PEdge> huang(log_threshold, log_sent_prob, total_out);

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

      pred_beam_huang_unary<typename Types::LBEdge> huang(log_threshold, log_sent_prob, total_out);

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
  auto c_r_s = [&](PEdge & edge, int i)->bool {
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
      edge.get_annotations() = std::move(a);
      return true;
    }
    return false;
  };
  
  for(size_t i=0; i<max_size; ++i) {
    Edge & edge = edges[i];
    if (c_r_s(edge. uedge(), i)) edge. uedge().replace_rule_probabilities(0);
    if (c_r_s(edge.lbedge(), i)) edge.lbedge().replace_rule_probabilities(0);
  }
}



template<class Types>
void PCKYAllCell<Types>::change_rules_resize(unsigned new_size, unsigned finer_idx)
{
  apply_on_lbedges(function<void(LBEdge&)>([new_size,finer_idx](LBEdge&e){
    //resize
    e.get_annotations().reset_probabilities(0.0);
    e.get_annotations().resize(new_size);
    //replace rule
    e.replace_rule_probabilities(finer_idx);
  })
  );
  apply_on_uedges(function<void(UEdge&)>([new_size,finer_idx](UEdge&e){
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
void PCKYAllCell<Types>::dump(std::ostream& out) const
{
  out << "(cell: span=" << get_end() - get_begin() << ", beg=" << get_begin() << " :"<< std::endl;
  int nb_entries = 0;
  for(unsigned i = 0; i < max_size ; ++i)
    if(not edges[i].is_closed()) {
      ++nb_entries;
      out << " " << i << ":" << edges[i] << std::endl;
    }
  out << "filled entries: " << nb_entries << ")";
}





#endif //PCKYALLCELL_HPP
