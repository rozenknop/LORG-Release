// -*- mode: c++ -*-
#ifndef _PARSERCKYALLMINDIVKB_CPP_
#define _PARSERCKYALLMINDIVKB_CPP_

#include "ParserCKYAllMinDivKB.h"
#include <utils/tick_count.h>
#include "utils/threads.h"




ParserCKYAllMinDivKB::ParserCKYAllMinDivKB(std::vector<AGrammar*>& cgs,
                                             const std::vector<double>& p, double b_t,
                                             const annot_descendants_type& annot_descendants_,
                                             bool accurate_, unsigned min_beam, int stubborn, unsigned k_)
: ParserCKYAll_Impl<MinDivKBTypes>(cgs, p, b_t, annot_descendants_, accurate_, min_beam, stubborn) , k(k_)
{
  // this is not in the super class because maxn parsing uses a
  //different mapping
  //create the coarse-to-fine map
  this->create_coarse_to_fine_mapping(this->grammars);

  UEdge::set_unary_chains(this->grammars[this->grammars.size() - 1]->get_unary_decoding_paths());
}


void ParserCKYAllMinDivKB::extract_solution()
{
  //  std::cout << "in extract" << std::endl;

  // this calls ParserCKYAllMinDivKB::compute_outside_probabilities()
  // and computes marginals for each daughter
  {
    //     BLOCKTIMING("compute_inside_outside_probabilities");
    compute_inside_outside_probabilities();
  }
//   // initialises q_insides and q_outsides to 1
//   function<void(MinDivProbabilityKB&)> reinit_1 = [](MinDivProbabilityKB&p){p.reinit_inside_outside(1.0);};
//   this->chart->opencells_apply([&](Cell & cell)
//   {
//     cell.apply_on_edges(reinit_1);
//   });
//   
//   // computes q from marginal probabilities on p
//   update_q();
  
  /* min divergence computation here ! : fixed point iterations */
  compute_inout_p();
  
  std::clog << *chart << std::endl;

  for (int i=0; i<1000; ++i) {
    //   while (false) {
      compute_inside_outside_q_probabilities();
//       std::clog << *chart << std::endl;
      compute_kl_distance();
      update_q();
  }
//   std::clog << *chart << std::endl;
  
  // compute the probability of the best subtree starting at each edge daughter
  // and initialize structures "candidates" and "derivations" of the edge
  {
    //     BLOCKTIMING("initialise_candidates");
    initialise_candidates();
  }
  //  std::cout << "after init cand" << std::endl;

  // complete structures "derivations" of the edge
  // from structures "candidates"
  {
    //     BLOCKTIMING("extend_all_derivations");
    extend_all_derivations();
  }
}


#include "parsers/ParserCKYAll.hpp"
#include "MinDivDaughters.hpp"

double MinDivProbabilityKB::log_normalisation_factor = 0;
double MinDivProbabilityKB::normalisation_factor = 0;
double MinDivProbabilityKB::log_normalisation_factor_q = 0;
double MinDivProbabilityKB::normalisation_factor_q = 0;
unsigned MinDivBest::size = 0;



ostream & MinDivProbabilityKB::operator>> (ostream & out) const 
{return out;}

inline void
MinDivBest::extend_derivation (unsigned int i, bool licence_unaries)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability_with_index& last = derivations[derivations.size() -1];

    //    std::cout << "last.probability " << last.probability << std::endl;

    assert(last.probability <= 0);

    find_succ(last,licence_unaries);
    //    std::cout << "after find_succ" << std::endl;
  }

  if (!candidates.empty()) {

    //get next element from the candidatesidates and append it to derivations
    pop_heap(candidates.begin(),candidates.end());
    derivations.push_back(candidates.back());
    candidates.pop_back();

    //    std::cout << "in edge " << edge << " there are " << derivations.size() << " derivations." << std::endl;

  }



#ifndef NDEBUG
   for(unsigned j = 0; j < derivations.size(); ++j) {
  //   std::cout << "derivations " << j << ": " << derivations[j].probability << " " ;
  //   if(!(derivations[j].dtrs)) { std::cout << "NULL"; }
  //   else {
  //     if(derivations[j].dtrs->is_lexical())
  //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //     if(derivations[j].dtrs->is_binary())
  //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //     if(!derivations[j].dtrs->is_binary() && !derivations[j].dtrs->is_lexical())
  //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(derivations[j].dtrs)->get_rule());
  //   }
  //   std::cout << std::endl;

    assert(derivations[j].probability <= 0);
  }
#endif

}

inline void MinDivBest::find_succ(packed_edge_probability_with_index& pep, bool licence_unaries)
{
  if(pep.dtrs->is_lexical())  { return;}
  // binary -> extend left and right daughters
  if(pep.dtrs->is_binary()) {
    const BinaryDaughter * d = static_cast<const BinaryDaughter*>(pep.dtrs);

    //extend to the left
    PEdge& left  = d->left_pdaughter() ;
    unsigned nextleft = pep.left_index + 1;
    left.get_best().extend_derivation(nextleft+1,true);

    // we haven't reached the expected number of solutions
    if(nextleft < left.get_best().n_deriv()) {

      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = d->tree_log_proba(p.left_index, p.right_index); 

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }

    }

    //extend to the right
    PEdge& right = d->right_pdaughter();
    unsigned nextright = pep.right_index + 1;

    right.get_best().extend_derivation(nextright+1,true);

    if(nextright < right.get_best().n_deriv()) {
      //        std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability_with_index p(pep);
      p.right_index = nextright;
      p.probability = d->tree_log_proba(p.left_index, p.right_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }

  //unary
  else {
    if(!licence_unaries) return;

    //      std::cout << "unary case" << std::endl;

    const UnaryDaughter* d = static_cast<const UnaryDaughter*>(pep.dtrs);

    //        std::cout << * d->get_rule() << std::endl;


    //extend to the left
    LBEdge& left  = d->lbdaughter();
    unsigned nextleft = pep.left_index + 1;

    left.get_best().extend_derivation(nextleft+1, false);

    if(nextleft < left.get_best().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = d->tree_log_proba(p.left_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;


      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }
}



inline void ParserCKYAllMinDivKB::extend_all_derivations()
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  Cell& root = chart->get_root();

  if (!root.exists_edge(start_symbol))
    //   //   std::cout << "no axiom at root" << std::endl;
    return;

  for (unsigned i = 2; i <= k; ++i)
  {
    //      std::cout << "before extend" << std::endl;
    root.get_edge(start_symbol).uedge().get_best().extend_derivation(i,true);
  }
}

/*****************************************************************************/
/*    real parsing stuff                                                     */
/*****************************************************************************/

ATOMIC_DOUBLE MinDivProbabilityKB::kl_distance_pq;

void ParserCKYAllMinDivKB::compute_kl_distance()
{
  MinDivProbabilityKB::kl_distance_pq = 0.0 ;
  chart->opencells_apply_top_down_nothread( [](Cell & cell){ cell.apply_on_edges(std::function<void(Edge&)>(MinDivProbabilityKB::update_kl_distance) ); } );
  MinDivProbabilityKB::kl_distance_pq += log2(MinDivProbabilityKB::normalisation_factor_q) /* * MinDivProbabilityKB::normalisation_factor*/;
  std::clog.precision(30); std::clog << "Divergence(p,q) = " << MinDivProbabilityKB::kl_distance_pq << std::endl;
}

void ParserCKYAllMinDivKB::compute_inout_p()
{
  chart->opencells_apply( [](Cell & cell){ cell.apply_on_edges(std::function<void(Edge&)>(MinDivProbabilityKB::compute_inout_p) ); } );
}

inline void ParserCKYAllMinDivKB::compute_outside_probabilities()
{
  MinDivProbabilityKB::set_normalisation_factor(get_sentence_probability());
  
  this->chart->opencells_apply_top_down([&](Cell & cell)
  {
      cell.apply_on_uedges (&       UnaryDaughter  ::outside_and_marginal);
      cell.apply_on_lbedges(&       BinaryDaughter ::outside_and_marginal, 
                            &       LexicalDaughter::outside_and_marginal);
  });
}

/****************************************************/
/*     computation of q insides                     */
/****************************************************/
inline void MinDivProbabilityKB::update_inside_q_lexical(const LexicalDaughter& dtr)
{
  inside_q += dtr.q ;
}
inline void MinDivProbabilityKB::update_inside_q_binary(const BinaryDaughter& dtr)
{
  inside_q += (
    dtr.q 
    * dtr.left_pdaughter().get_prob_model().inside_q
    * dtr.right_pdaughter().get_prob_model().inside_q
  );
}
inline void MinDivProbabilityKB::update_inside_q_unary(const UnaryDaughter& dtr)
{
  inside_q += (
    dtr.q
    * dtr.lbdaughter().get_prob_model().inside_q
    );
}

inline void ParserCKYAllMinDivKB::compute_inside_q_probabilities()
{
  this->chart->opencells_apply_bottom_up([&](Cell & cell)
  {
    cell.apply_on_lbedges(& MinDivProbabilityKB::update_inside_q_lexical,
                          & MinDivProbabilityKB::update_inside_q_binary);
    
    cell.apply_on_uedges(& MinDivProbabilityKB::update_inside_q_unary);
  }
  );

  MinDivProbabilityKB::set_normalisation_factor_q(get_sentence_probability_q());
}

/****************************************************/
/*     computation of q outsides                    */
/****************************************************/


inline void MinDivProbabilityKB::update_outside_q_binary(BinaryDaughter& dtr)
{
  #ifdef USE_THREADS
  tbb::atomic<double> * left  = (tbb::atomic<double> *) &dtr. left_pdaughter().get_prob_model().outside_q;
  tbb::atomic<double> * right = (tbb::atomic<double> *) &dtr.right_pdaughter().get_prob_model().outside_q;
  #else
  double * left  = &dtr. left_pdaughter().get_prob_model().outside_prob;
  double * right = &dtr.right_pdaughter().get_prob_model().outside_prob;
  #endif
  *left += (
    outside_q
    * dtr.q 
    * dtr.right_pdaughter().get_prob_model().inside_q
  );
  *right += (
    outside_q
    * dtr.q 
    * dtr.left_pdaughter().get_prob_model().inside_q
  );
  dtr.mq = dtr.q * outside_q *  dtr.left_pdaughter().get_prob_model().inside_q * dtr.right_pdaughter().get_prob_model().inside_q / get_normalisation_factor_q();
}
inline void MinDivProbabilityKB::update_outside_q_unary(UnaryDaughter& dtr)
{
  dtr. lbdaughter().get_prob_model().outside_q += (
    outside_q
    * dtr.q
  );
  dtr.mq = dtr.q * outside_q * dtr. lbdaughter().get_prob_model().inside_q  / get_normalisation_factor_q();
}

inline void MinDivProbabilityKB::update_outside_q_lexical(LexicalDaughter& dtr)
{
  dtr.mq = dtr.q * outside_q / get_normalisation_factor_q();
}


inline void ParserCKYAllMinDivKB::compute_outside_q_probabilities()
{
  this->chart->opencells_apply_top_down([&](Cell & cell)
  {
    cell.apply_on_uedges(& MinDivProbabilityKB::update_outside_q_unary);
    cell.apply_on_lbedges(& MinDivProbabilityKB::update_outside_q_binary);
    cell.apply_on_lbedges(& MinDivProbabilityKB::update_outside_q_lexical);
  }
  );
}


/*********************************************************************/
/*    insides and outsides on q                                      */
/*********************************************************************/
inline void MinDivProbabilityKB::reinit_inside_outside(double val)
{
  inside_q = outside_q = inside_p = outside_p = val;
}

inline void MinDivProbabilityKB::reinit_q_inside_outside(double val)
{
  inside_q = outside_q = val;
}

inline void ParserCKYAllMinDivKB::compute_inside_outside_q_probabilities()
{
  function<void(MinDivProbabilityKB&)> reinit_0 = [](MinDivProbabilityKB&p){p.reinit_q_inside_outside(0);};
  this->chart->opencells_apply([&](Cell & cell)
  {
    cell.apply_on_edges(reinit_0);
  });
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);
  this->chart->get_root().get_edge(start_symbol).uedge().get_prob_model().set_outside_q(1.0);
  
  compute_inside_q_probabilities();

  compute_outside_q_probabilities();
}

double ParserCKYAllMinDivKB::get_sentence_probability_q() const
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  if(chart->get_root().exists_uedge(start_symbol))
    return chart->get_root().get_edge(start_symbol).uedge().get_prob_model().get_inside_q();
  else
    return LorgConstants::NullProba;
}

/***********************************************************************/
/*      update q as marginal(p) * sentence_probability(q) / (inside(q)*outside(q))               */
/***********************************************************************/

inline double formula(double q, double mp, double sq, double ioq, double iop)
{
  return mp / iop ;
//   return mp / ioq ;
//     return std::min(1.,mp / ioq) ;
//     return mp * sq / ioq;
//    return std::min(1.0, mp * sq / ioq) ;
//     return mp<1 ? (mp/(1-mp))*(sq/ioq - q) : (sq/ioq);
//     return std::min(1., mp<1 ? (mp/(1-mp))*(sq/ioq - q) : (sq/ioq));
}

inline void MinDivProbabilityKB::update_q_lexical(LexicalDaughter& dtr)
{
  if (outside_q != LorgConstants::NullProba)
//     dtr.q = dtr.mp / outside_prob;
//     dtr.q += (dtr.mp / outside_prob - dtr.q) * 0.01;
//      dtr.q += (dtr.mp / outside_q - dtr.q) * 0.01;
    dtr.q = formula(dtr.q, dtr.mp, get_normalisation_factor_q(), outside_q, outside_p/normalisation_factor) ;
  else
    dtr.q = 0 ;
}

inline void MinDivProbabilityKB::update_q_unary(UnaryDaughter& dtr)
{
  if (outside_q != LorgConstants::NullProba)
    dtr.q = formula(dtr.q, dtr.mp, get_normalisation_factor_q(),
                    outside_q * dtr.lbdaughter().get_prob_model().inside_q,
                    outside_p * dtr.lbdaughter().get_prob_model().inside_p / normalisation_factor
                   ) ;
//     dtr.q = dtr.mp * get_normalisation_factor_q() / (
//       outside_q
//       * dtr.lbdaughter().get_prob_model().inside_q
//     );
//     dtr.q = dtr.mp / (
//       outside_prob
//       * dtr.lbdaughter().get_prob_model().inside_prob
//     );
//     dtr.q += (dtr.mp / (outside_q * dtr.lbdaughter().get_prob_model().inside_q) - dtr.q) * 0.01;
}
inline void MinDivProbabilityKB::update_q_binary(BinaryDaughter& dtr)
{
  if (outside_q != LorgConstants::NullProba)
    dtr.q = formula(dtr.q, dtr.mp, get_normalisation_factor_q(),
                    outside_q * dtr. left_pdaughter().get_prob_model().inside_q * dtr.right_pdaughter().get_prob_model().inside_q,
                    outside_p * dtr. left_pdaughter().get_prob_model().inside_p * dtr.right_pdaughter().get_prob_model().inside_p / normalisation_factor) ;
//     dtr.q = dtr.mp * get_normalisation_factor_q() / (
//       outside_q
//       * dtr. left_pdaughter().get_prob_model().inside_q
//       * dtr.right_pdaughter().get_prob_model().inside_q
//     );
//     dtr.q = dtr.mp / (
//       outside_prob
//       * dtr. left_pdaughter().get_prob_model().inside_prob
//       * dtr.right_pdaughter().get_prob_model().inside_prob
//     );
//     dtr.q += (dtr.mp / (outside_q * dtr. left_pdaughter().get_prob_model().inside_q * dtr.right_pdaughter().get_prob_model().inside_q) - dtr.q) * 0.01;
}


inline void ParserCKYAllMinDivKB::update_q()
{

//   this->chart->opencells_apply([&](Cell & cell)
//   {
//     cell.apply_on_lbedges(& MinDivProbabilityKB::update_q_lexical,
//                           & MinDivProbabilityKB::update_q_binary);
//     cell.apply_on_uedges(& MinDivProbabilityKB::update_q_unary);
//   });
  this->chart->opencells_apply_bottom_up([&](Cell & cell)
  {
    cell.apply_on_edges( function<void(Edge &)>([&](Edge & e){
      if (e.lbedge().is_opened()) {
        for (auto & dtr: e.lbedge().get_lexical_daughters()) {
          e.lbedge().get_prob_model().update_q_lexical(dtr);
          compute_inside_outside_q_probabilities();
        }
        for (auto & dtr: e.lbedge().get_binary_daughters()) {
          e.lbedge().get_prob_model().update_q_binary(dtr);
          compute_inside_outside_q_probabilities();
        }
      }
      if (e.uedge().is_opened()) {
        for (auto & dtr: e.uedge().get_unary_daughters()) {
          e.uedge().get_prob_model().update_q_unary(dtr);
          compute_inside_outside_q_probabilities();
        }
      }
    }));
  });
}


/**************************************************************/
/* Computation of the Kullback-Leibler distance KL(p,q)       */
/**************************************************************/

inline void MinDivProbabilityKB::update_kl_distance(Edge & edge)
{
  if (edge.uedge().is_opened()) {
    for (const auto & dtr : edge.uedge().get_unary_daughters()) {
      kl_distance_pq += dtr.get_rule()->entropy_term(edge.get_annotations().outside_probabilities.array,
                                                     dtr.lbdaughter().get_annotations().inside_probabilities.array) / normalisation_factor
      - dtr.mp * log2(dtr.q);
    }
  }
  if (edge.lbedge().is_opened()) {
    for (const auto & dtr : edge.lbedge().get_lexical_daughters()) {
      kl_distance_pq += dtr.get_rule()->entropy_term(edge.get_annotations().outside_probabilities.array) / normalisation_factor
      - dtr.mp * log2(dtr.q);
    }
    for (const auto & dtr : edge.lbedge().get_binary_daughters()) {
      kl_distance_pq += dtr.get_rule()->entropy_term(edge.get_annotations().outside_probabilities.array,
                                                     dtr.left_pdaughter().get_annotations().inside_probabilities.array,
                                                     dtr.right_pdaughter().get_annotations().inside_probabilities.array) / normalisation_factor
      - dtr.mp * log2(dtr.q);
    }
  }
}

void MinDivProbabilityKB::compute_inout_p(MinDivProbabilityKB::Edge& edge)
{
  if (edge.uedge().is_opened()) {
    auto & e = edge.uedge();
    auto & pr = e.get_prob_model();
    pr.inside_p = pr.outside_p = 0;
    for (auto & p : e.get_annotations().inside_probabilities.array) {
      if (p!=LorgConstants::NullProba) pr.inside_p += p;
    }
    for (auto & p : e.get_annotations().outside_probabilities.array) {
      if (p!=LorgConstants::NullProba) pr.outside_p += p;
    }
  }
  if (edge.lbedge().is_opened()) {
    auto & e = edge.lbedge();
    auto & pr = e.get_prob_model();
    pr.inside_p = pr.outside_p = 0;
    for (auto & p : e.get_annotations().inside_probabilities.array) {
      if (p!=LorgConstants::NullProba) pr.inside_p += p;
    }
    for (auto & p : e.get_annotations().outside_probabilities.array) {
      if (p!=LorgConstants::NullProba) pr.outside_p += p;
    }
  }
}

/**************************************************************/
/* Filling the structures before Best Tree extraction         */
/**************************************************************/

template<class TDaughter>
inline void MinDivBest::update_best(const TDaughter& dtr)
{
//   std::cout << "updating "<<*this << "(dtr: " << *dtr.get_rule() << ")" << std::endl;
  packed_edge_probability_with_index pep;

  pep.probability = dtr.tree_log_proba();

//   std::cout << "proba= " << pep.probability << std::endl;
  assert(pep.probability <=0);

  pep.dtrs = &dtr;

  candidates.push_back(pep);

  if(derivations.empty())
    derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
    derivations[0] = pep;

//   std::cout << "updated  "<<*this << std::endl;
}

inline void MinDivBest:: finalize_best()
{
  if(!candidates.empty()) {
    if(candidates.size() > size) {
        std::nth_element(candidates.begin(),candidates.begin()+size,candidates.end(),std::greater<packed_edge_probability>());
        candidates.resize(size);
  }
  std::make_heap(candidates.begin(),candidates.end());
  
  std::pop_heap(candidates.begin(),candidates.end());
  candidates.pop_back();
}
}

#include <ostream>

inline void ParserCKYAllMinDivKB::fill_bests()
{
//   std::cout << * this->chart << std::endl;

  this->chart->opencells_apply_bottom_up(
    [](Cell & cell)
    {
//       std::cout << "filling cell " << &cell << " : ======================================================" << cell << std::endl;
      cell.apply_on_lbedges (/*function<void(Edge&)>([](Edge&e){(std::cout << "(edge.1:"<<&e <<" : " << e << ") ").flush();}), */
                             toFunc(&Best::update_best<LexicalDaughter>),
//                            function<void(Edge&)>([](Edge&e){(std::cout << "(edge.2:"<<&e <<" : " << e << ") ").flush();}), 
                             toFunc (&Best::update_best<BinaryDaughter>));
      cell.apply_on_uedges (toFunc(&Best::update_best<UnaryDaughter>));
      cell.apply_on_edges (toFunc (&Best::finalize_best));
//       std::cout << "best filled for cell " << &cell << " : " << cell << std::endl;
    }
  );
}


inline void ParserCKYAllMinDivKB::initialise_candidates()
{
  
  double sentence_probability = std::log(get_sentence_probability());
  //  unsigned sent_size = chart->get_size();
  
  MinDivProbabilityKB::set_log_normalisation_factor(sentence_probability);
  MinDivBest::set_size(k);
  
  fill_bests();
}




template<>
std::ostream & operator<<(std::ostream & out, const BasePackedEdge<MinDivKBTypes> & e)
{
  out << "  (in,out)=(" << e.get_prob_model().get_inside_q() << ", " << e.get_prob_model().get_outside_q() << ")";
  return out;
}

template<>
std::ostream & operator<<(std::ostream & out, const UPackedEdge<MinDivKBTypes> & e)
{
  out << (BasePackedEdge<MinDivKBTypes>&)e;
  if (not e.get_unary_daughters().empty()) {
    out << std::endl << "  unary dtrs : ";
    for (const auto& dtr: e.get_unary_daughters()) {
      out << " " << SymbolTable::instance_nt().translate(dtr.get_rule()->get_lhs());
      out << "->" << SymbolTable::instance_nt().translate(dtr.get_rule()->get_rhs0());
      out << " (" << dtr.q << " ["<< dtr.mq << " ~ " << dtr.mp << "])";
    }
  }
  return out;
}

template<>
std::ostream & operator<<(std::ostream & out, const LBPackedEdge<MinDivKBTypes> & e)
{
  out << (BasePackedEdge<MinDivKBTypes>&)e;
  if (not e.get_lexical_daughters().empty()) {
    out << std::endl << "  lexical dtrs : ";
    for (const auto & dtr: e.get_lexical_daughters()) {
      out << " " << SymbolTable::instance_nt().translate(dtr.get_rule()->get_lhs());
      out << "->" << SymbolTable::instance_word().translate(dtr.get_rule()->get_rhs0());
      out << " (" << dtr.q << " ["<< dtr.mq << " ~ " << dtr.mp << "])";
    }
  }
  if (not e.get_binary_daughters().empty()) {
    out << std::endl << "  binary dtrs : ";
    for (const auto & dtr: e.get_binary_daughters()) {
      out << " " << SymbolTable::instance_nt().translate(dtr.get_rule()->get_lhs());
      out << "->" << SymbolTable::instance_nt().translate(dtr.get_rule()->get_rhs0());
      out << " " << SymbolTable::instance_nt().translate(dtr.get_rule()->get_rhs1());
      out << " (" << dtr.q << " ["<< dtr.mq << " ~ " << dtr.mp << "])";
    }
  }
  return out;
}

template<>
std::ostream & operator<<(std::ostream & out, const PackedEdge<MinDivKBTypes> & e)
{
  out << "(edge: " << &e << " " << (AnnotatedEdge<MinDivKBTypes>&)(e);
  if (e.uedge().is_opened()) out << std::endl << e.uedge() ;
  if (e.lbedge().is_opened()) out << std::endl << e.lbedge() ;
  out << ")";
  return out;
}




#endif /* _PARSERCKYALLMINDIVKB_CPP_ */
