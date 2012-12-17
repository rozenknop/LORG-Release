#ifndef ParserCKYAll_hpp
#define ParserCKYAll_hpp

#include "ParserCKYAll.h"
#include "ChartCKY.hpp"

template <class Types>
ParserCKYAll_Impl<Types>::ParserCKYAll_Impl(std::vector<AGrammar*>& cgs,
                                            const std::vector<double>& p,
                                            double prior_threshold,
                                            const annot_descendants_type& annot_descendants_,
                                            bool accurate_,
                                            unsigned min_beam, int stubborn) :
    ParserCKYAll(cgs, p, prior_threshold, annot_descendants_, accurate_, min_beam, stubborn),
  chart(NULL)
{};



template <class Types>
ParserCKYAll_Impl<Types>::~ParserCKYAll_Impl()
{
  for (std::vector<AGrammar*>::iterator i(grammars.begin()); i != grammars.end(); ++i)
    if(i != grammars.begin()) // the first grammar is deleted by super class
    {
      delete *i;
      *i = NULL;
    }
}

#include "utils/tick_count.h"

template <class Types>
void ParserCKYAll_Impl<Types>::parse(int start_symbol) const
{
  int ntries = stubbornness;
  double beam_threshold = prior_beam_threshold;

  do {
    
    //clear only when first try was a failure
    if(ntries != stubbornness) {
      chart->prepare_retry();
    }

    // last resort
    if(ntries == 0)
      beam_threshold = 0;

//        std::clog << "ParserCKY::parse ntries = " << ntries << " threshold : " << beam_threshold << std::endl;


    //init
    {
      //               BLOCKTIMING("parse_init");
      bool beam_short = chart->get_size() >= min_length_beam;
      chart->opencells_apply([&](Cell& cell){
        if(!cell.is_empty()) {
          this->add_unary_init(cell,cell.get_top());
          //           std::cout << cell << std::endl;
          cell.apply_on_uedges( &Edge::add_unary_insides_to_merged );

          // prevent short sentences from being skipped ...
          if(beam_short)
            cell.beam(priors, beam_threshold);

          // if(cell.is_closed())
            //   std::cout << "(" << i << "," <<j << ") is closed" << std::endl;
//           std::cout << cell << std::endl;
        }
      }
      );
    }
//     std::clog << "ParserCKY::parse after init" << std::endl;

    //actual cky is here
    {
      //               BLOCKTIMING("process_internal_rules");
      process_internal_rules(beam_threshold);
    }
//     std::clog << "ParserCKY::parse after process_internal_rules(beam_threshold)" << std::endl;

    if(ntries == 0)
      break;

    --ntries;
    beam_threshold /= 10;
  }
  while (stubbornness >=0 &&
         beam_threshold > 0 &&
         !chart->get_root().exists_uedge(start_symbol));
}

template <class Types>
void ParserCKYAll_Impl<Types>::process_internal_rules(double beam_threshold) const
{
  chart->opencells_apply_bottom_up(
    [&,beam_threshold](Cell&cell)
    {
//       std::clog << "ParserCKYAll_Impl<Types>::process_internal_rules : before process_cell(" << cell.get_end() - cell.get_begin() + 1 << ","<< cell.get_begin() << ") (" << &cell << ")" << std::endl;
//       std::clog << cell << std::endl;
      this->process_cell(cell, beam_threshold);
    },
    1 // start from span = 1 (i.e. 2 words !)
  );
}

template <class Types>
void ParserCKYAll_Impl<Types>::process_cell(Cell& cell, double beam_threshold) const
{
//   BLOCKTIMING("ParserCKYAll_Impl<Types>::process_cell");
  const unsigned & begin = cell.get_begin();
  const unsigned & end   = cell.get_end();
  const bool & isroot = cell.get_top();

  // look for all possible new edges

  //application of binary rules
  {
    // BLOCKTIMING("process_cell binary");
    for (unsigned m = begin; m < end; ++m) {
      // m is the mid-point
      Cell& left_cell = chart->access(begin,m);
      if(!left_cell.is_closed()) {
        Cell& right_cell = chart->access(m+1,end);
        if( !right_cell.is_closed()) {
          get_candidates(left_cell,right_cell,cell);
          
        }
      }
    }
//     std::clog << "ParserCKYAll_Impl<Types>::process_cell : " << cell << std::endl;
  }
  //unary rules
  {
    // BLOCKTIMING("process_cell unary");
    add_unary_internal(cell, isroot);
  }
  cell.apply_on_uedges( &Edge::add_unary_insides_to_merged );

  //   std::clog << "ParserCKYAll_Impl<Types>::process_cell : after add_unary_internal : " << cell << std::endl;
  // pruning
  if(chart->get_size() >= min_length_beam)
  {
    // BLOCKTIMING("process_cell beam");
    cell.beam(priors, beam_threshold);
//     std::clog << "ParserCKYAll_Impl<Types>::process_cell : after cell.beam : " << cell << std::endl;
  }
  // if(cell.is_closed())
  //   std::cout << "(" << begin << "," << end << ") is closed" << std::endl;
}

template <class Types>
inline
void ParserCKYAll_Impl<Types>::get_candidates(Cell& left_cell,
                                              Cell& right_cell,
                                              Cell& result_cell) const
{
  BLOCKTIMING("ParserCKYAll_Impl<Types>::get_candidates (binaires)");
//   std::clog << "ParserCKYAll_Impl<Types>::get_candidates("<<&left_cell<<","<<&right_cell<<","<<&result_cell<<")"<<std::endl;
  //   {
    //               BLOCKTIMING("get_candidates counting");
    // count the number of daughters to create
    //     std::vector<int> nb_rules(result_cell.get_max_size(), 0);
    //     for (const auto & same_rhs0_rules: brules) {
      //       if (left_cell.exists_edge(same_rhs0_rules.rhs0)) {
        //         for(const auto & same_rhs: same_rhs0_rules) {
          //           if (right_cell.exists_edge(same_rhs.rhs1)) {
            //             for(const auto & rule: same_rhs) {
              //               ++ nb_rules[rule->get_lhs()];
            //             }
            //           }
            //         }
            //       }
            //     }
            //     // create daughters
            //     result_cell.reserve_binary_daughters(nb_rules);
            //   }  
  {
    //               BLOCKTIMING("get_candidates creating");
    //iterating through all the rules P -> L R, indexed by L
    for (const auto & same_rhs0_rules: brules) {
      Edge & left_edge = left_cell.get_edge(same_rhs0_rules.rhs0) ;
      if (left_edge.is_opened()) {
        const double & L = left_edge.get_annotations().inside_probabilities.array[0];
        //iterating through all the rules P -> L R, indexed by R, L fixed
        for(const auto & same_rhs: same_rhs0_rules) {
          Edge & right_edge = right_cell.get_edge(same_rhs.rhs1);
          if (right_edge.is_opened()) {
            double LR = L * right_edge.get_annotations().inside_probabilities.array[0];
            
            //iterating through all the rules P -> L R, indexed by P, R and L fixed
            for(const auto & rule: same_rhs) {
              BLOCKTIMING("ParserCKYAll_Impl<Types>::get_candidates result_cell.process_candidate(*left_edge,*right_edge, rule, LR);");
//                   std::clog << "calling result_cell.process_candidate(" << left_edge << " , " << right_edge << " , " << rule << " , " << LR << ")" << std::endl;
              result_cell.process_candidate(left_edge,right_edge, rule, LR);
            }
          }
        }
      }
    }
  }
}


template <class Types>
inline
void ParserCKYAll_Impl<Types>::add_unary_init(Cell& cell, bool isroot) const
{
  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs]
  static std::vector<short>::const_iterator unary_rhs_itr_begin = unary_rhs_from_pos.begin();
  static std::vector<short>::const_iterator unary_rhs_itr_end = unary_rhs_from_pos.end();

  for(std::vector<short>::const_iterator unary_rhs_itr(unary_rhs_itr_begin); unary_rhs_itr != unary_rhs_itr_end; ++unary_rhs_itr) {

    if (cell.exists_lbedge(*unary_rhs_itr)) {
//       BLOCKTIMING("ParserCKYAll_Impl<Types>::add_unary_init");
      process_unary(cell,*unary_rhs_itr, isroot);
    }
  }
}

template <class Types>
inline
void ParserCKYAll_Impl<Types>::add_unary_internal(Cell& cell, bool isroot) const
{
  //for each unary rule set in the grammar [sets made up of all unary rules with a particular rhs being a lhs of a binary rule]
  std::vector<short>::const_iterator unary_rhs_itr_end = unary_rhs_from_binary.end();
  for(std::vector<short>::const_iterator unary_rhs_itr = unary_rhs_from_binary.begin();unary_rhs_itr!=unary_rhs_itr_end;++unary_rhs_itr) {

    if (cell.exists_lbedge(*unary_rhs_itr)) {
      //BLOCKTIMING("ParserCKYAll_Impl<Types>::add_unary_internal");
      process_unary(cell,*unary_rhs_itr,isroot);
    }
  }
}


template <class Cell>
struct processunary
{
  Cell& cell;
  double R_inside;
  processunary(Cell& c, double L) : cell(c), R_inside(L) {};
  void operator()(const URuleC2f* r) const
  {
    cell.process_candidate((typename Cell::UnaryRule *)r,R_inside);
  }
};


template <class Types>
void ParserCKYAll_Impl<Types>::process_unary(Cell& cell, int rhs, bool isroot) const
{
  //BLOCKTIMING("ParserCKYAll_Impl<Types>::process_unary");
  const std::vector<const URuleC2f*>& rules = isroot ?
                                              unary_rhs_2_rules_toponly[rhs] :
                                              unary_rhs_2_rules_notop[rhs];

  double R_inside = cell.get_edge(rhs).get_annotations().inside_probabilities.array[0];

  std::for_each(rules.begin(),rules.end(),processunary<Cell>(cell, R_inside));
}



template <class Types>
void ParserCKYAll_Impl<Types>::compute_merged_outside_probabilities()
{
  this->chart->opencells_apply_top_down( & Cell::compute_merged_outside_probabilities) ;
}

template <class Types>
void ParserCKYAll_Impl<Types>::compute_merged_inside_probabilities()
{
  this->chart->opencells_apply_bottom_up( & Cell::compute_merged_inside_probabilities );
}


template <class Types>
void ParserCKYAll_Impl<Types>::edge_daughters_differentiation()
{
  this->chart->opencells_apply(std::function<void(Cell&)>([](Cell&cell){
      cell.apply_on_edges(&Edge::daughters_differentiation);
  }));
}

template <class Types>
void ParserCKYAll_Impl<Types>::compute_outside_probabilities()
{
  this->chart->opencells_apply_top_down( & Cell::compute_outside_probabilities) ;
}

template <class Types>
void ParserCKYAll_Impl<Types>::compute_inside_probabilities()
{
  this->chart->opencells_apply_bottom_up( & Cell::compute_inside_probabilities );
}


template <class Types>
double ParserCKYAll_Impl<Types>::get_sentence_probability() const
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);

  if(chart->get_root().exists_uedge(start_symbol))
    return chart->get_root().get_edge(start_symbol).uedge().get_annotations().get_inside(0);
  else
    return LorgConstants::NullProba;
}


//absolute beam
template <class Types>
void ParserCKYAll_Impl<Types>::beam_chart(double log_sent_prob, double log_threshold, bool huang)
{
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);
//   std::clog << *chart << std::endl;
  chart->get_root().get_edge(start_symbol).get_annotations().reset_outside_probabilities(1.0);
  compute_merged_outside_probabilities();
//   std::clog << "After compute_outside_probabilities :" << *chart << std::endl;
//   if (chart->get_root().is_closed() || !chart->get_root().exists_uedge(start_symbol))
//     std::cout << "no axiom at root after beam_chart.compute_outside_probabilities" << std::endl;

  this->chart->opencells_apply_bottom_up(
      [log_sent_prob, log_threshold, huang, &chart, &start_symbol]
      (Cell& cell)
      {
        cell.apply_on_lbedges(&LBEdge::clean_invalidated_binaries,
                              std::function<void(Edge&)>([](Edge&e){if (e.lbedge().no_daughters()) e.close_lb();}));
//         std::clog << "after clean_invalidated_binaries :" << *chart << std::endl;
//         if (chart->get_root().is_closed() || !chart->get_root().exists_uedge(start_symbol))
//           std::cout << "no axiom at root after beam_chart.close_lb(1) of cell (" << cell.get_end()-cell.get_begin()+1 << ","<<cell.get_begin()<<")"<<std::endl;
        cell.beam(log_threshold, log_sent_prob);
//         std::clog << "after beam :" << *chart << std::endl;
//         if (chart->get_root().is_closed() || !chart->get_root().exists_uedge(start_symbol))
//           std::cout << "no axiom at root after beam_chart.beam of cell (" << cell.get_end()-cell.get_begin()+1 << ","<<cell.get_begin()<<")"<<std::endl;
        cell.clean();
//         std::clog << "after clean :" << *chart << std::endl;
//         if (chart->get_root().is_closed() || !chart->get_root().exists_uedge(start_symbol))
//           std::cout << "no axiom at root after beam_chart.clean(1) of cell (" << cell.get_end()-cell.get_begin()+1 << ","<<cell.get_begin()<<")"<<std::endl;
        
        if(!cell.is_closed() && huang) {
          cell.apply_on_lbedges(&LBEdge::clean_invalidated_binaries,
                                std::function<void(Edge&)>([](Edge&e){if (e.lbedge().no_daughters()) e.close_lb();}));
          cell.beam_huang(std::log(0.0001), log_sent_prob);
          cell.clean();
        }
      }
  );
}




/////////////////////////////
//// mapping c2f ////////////
/////////////////////////////
// should be moved somewhere else


//  calculates c2f mapping
// returns rules that don't belong to the mapping
template <typename Key, typename MyRule>
std::vector<MyRule*> calculate_mapping(typename rulevect2mapvect<Key,MyRule>::map_type& map, unsigned size)
{
  std::vector<MyRule*> r2remove;
  for(typename rulevect2mapvect<Key,MyRule>::map_type::const_iterator i(map.begin()); i != map.end(); ++i)
  {
    if(i->second.size() != size
       || (std::find_if(i->second.begin(), i->second.end(), std::mem_fun(&MyRule::is_empty)) != i->second.end())
       )
      r2remove.push_back(i->second[0]);
    else
      for(unsigned g = 0 ; g < size - 1; ++g)
        i->second[g]->add_finer(i->second[g+1]);
  }
  return r2remove;
}



// calls previous function
// and removes useless rules
template <typename Key, typename MyRule>
void process_internal(typename rulevect2mapvect<Key,MyRule>::map_type& map, std::vector<MyRule>& grammar_coarse_rules, unsigned size)
{
  std::vector<MyRule*> r2remove = calculate_mapping<Key,MyRule>(map,size);
  typename std::vector<MyRule>::iterator end = grammar_coarse_rules.end();
  for(typename std::vector<MyRule*>::iterator i(r2remove.begin()); i != r2remove.end(); ++i) {
    // ++r
    // if(r % 1000 == 0)
    //   std::cout << "removing " << r << " of " << coarse_rules_wo_finers.size() << std::endl;
    end = std::remove(grammar_coarse_rules.begin(), end,**i);
  }
  grammar_coarse_rules.erase(end,grammar_coarse_rules.end());
}


#define MAP std::unordered_map
//#define MAP std::map

template <class Types>
void ParserCKYAll_Impl<Types>::create_coarse_to_fine_mapping(std::vector<AGrammar*>& cgs)
{
  //  std::clog << "before mapping" << std::endl;

  typedef std::pair<int, std::pair<int,int> > bkey;
  typedef std::pair<int,int> ukey;

  MAP< bkey, std::vector<BRuleC2f*> > bmap;
  MAP< ukey, std::vector<URuleC2f*> > umap;
  MAP< ukey, std::vector<LexicalRuleC2f*> > lmap;

  rulevect2mapvect<bkey,BRuleC2f> bc2f(bmap);
  rulevect2mapvect<ukey, URuleC2f> uc2f(umap);
  rulevect2mapvect<ukey, LexicalRuleC2f> lc2f(lmap);

  for(std::vector<AGrammar*>::const_iterator g(cgs.begin()); g != cgs.end(); ++g) {
    bc2f.add_all((*g)->binary_rules);
    uc2f.add_all((*g)->unary_rules);
    lc2f.add_all((*g)->lexical_rules);
  }

  process_internal<bkey,BRuleC2f>(bmap, cgs[0]->binary_rules, cgs.size());
  process_internal<ukey,URuleC2f>(umap, cgs[0]->unary_rules, cgs.size());

  std::vector<LexicalRuleC2f*> l2remove = calculate_mapping<ukey,LexicalRuleC2f>(lmap, cgs.size());
  for(std::vector<LexicalRuleC2f*>::iterator i(l2remove.begin()); i != l2remove.end(); ++i) {
    remove_lex_rule(*i);
  }

  //  std::clog << "after mapping" << std::endl;

}

////////////////////////////////
/////////////// C2f ///////////
///////////////////////////////
template <class Types>
void ParserCKYAll_Impl<Types>::beam_c2f(int start_symbol)
{
  if(chart->get_root().exists_uedge(start_symbol)) {
    beam_c2f(grammars, annot_descendants);
  }
}

template <class Types>
void ParserCKYAll_Impl<Types>::beam_c2f(const std::vector<AGrammar*>& current_grammars,
                                        const annot_descendants_type& /*current_annot_descendants*/)
{
  static int top_idx = SymbolTable::instance_nt().get_label_id(LorgConstants::tree_root_name);

  //  std::cout << "beam_c2f" << std::endl;

  for(unsigned i = 0; i < current_grammars.size() - 1; ++i) {

    double beam_threshold = io_beam_thresholds[i + 1];

    // std::cout << std::log(get_sentence_probability()) << std::endl;
    //     std::cout << "beaming with grammar: " << i << std::endl;


    // FIX: This test messes with product grammar parsing
    // TODO: Do this test only with the first grammar
//     if(i != 0) {// inside_probs already computed when bulding the chart
         std::clog << "ParserCKYAll_Impl<Types>::beam_c2f before inside, sentence probability = " << std::log(get_sentence_probability()) << std::endl;
         compute_merged_inside_probabilities();
//     }



    // if(chart->get_root().is_closed())
    //   std::cout << "root cell is closed" << std::endl;
    // else if(!chart->get_root().exists_edge(top_idx))
    //   std::cout << "top is not in root cell" << std::endl;

    if(chart->get_root().is_closed() || !chart->get_root().exists_uedge(top_idx)) {
//             std::cerr << "grammar " << i << " spoiled the fun :(" << std::endl;
      break;
    }
    //    std::cout << "after inside" << std::endl;
    //    std::cout << "before beam" << std::endl;
    double sp = std::log(get_sentence_probability());
    std::clog << "ParserCKYAll_Impl<Types>::beam_c2f after inside, sentence probability: " << sp << std::endl;

    // huang beam seems to affect only the first pass
    //bool huang = i == 0;
    bool huang = false;
    if(chart->get_size() >= min_length_beam) // TODO if sentence is short skip everything but correct resizing
      beam_chart(sp, beam_threshold, huang);
    //        std::cout << "after beam" << std::endl;

    // PCKYAllCell& root = chart->get_root();
    if (chart->get_root().is_closed() || !chart->get_root().exists_uedge(top_idx))
      std::cout << "no axiom at root after beam" << std::endl;


    //    std::cout << "before change" << std::endl;

    // TODO this function should take current_annot_descendants as an argument
    // instead annot_descendants is changed in ParserCKYAllMaxVarMultiple::extract_solution
    // which is a bit .. hackish
    change_rules_resize(i, current_grammars);
    //     std::clog << "after change_rules_resize :" << *chart << std::endl;
  }
}

template <class Types>
void ParserCKYAll_Impl<Types>::change_rules_resize(unsigned step,
                                                   const std::vector<AGrammar*>& current_grammars) const
{
  const AnnotatedLabelsInfo& next_annotations = current_grammars[step+1]->get_annotations_info();
  const std::vector<std::vector<std::vector<unsigned> > >& annot_descendants_current =  annot_descendants[step];

  this->chart->opencells_apply(
    [next_annotations, annot_descendants_current]
    (Cell& cell)
    {
      cell.change_rules_resize(next_annotations, annot_descendants_current);
    }
  );
}

template <class Types>
void ParserCKYAll_Impl<Types>::get_parses(int start_symbol, unsigned kbest,
                                          bool always_output_forms, bool output_annotations,
                                          std::vector<std::pair<PtbPsTree *,double> >& best_trees,
                                          bool unary_start_edge)
{
  for(unsigned i = 0; i < kbest; ++i) {
    // get results
    if(!chart->has_solution(start_symbol, i, unary_start_edge)) {
      break;
    }
    PtbPsTree * t = chart->get_best_tree(start_symbol, i, always_output_forms, output_annotations, unary_start_edge);
    best_trees.push_back(std::make_pair(t, chart->get_score(start_symbol, i, unary_start_edge)));
  }

}

template<class Types>
inline
typename ParserCKYAll_Impl<Types>::AGrammar& ParserCKYAll_Impl<Types>::get_grammar(unsigned idx)
{
  return *(grammars[idx]);
}

template<class Types>
inline
const typename ParserCKYAll_Impl<Types>::AGrammar& ParserCKYAll_Impl<Types>::get_grammar(unsigned idx) const
{
  return *(grammars[idx]);
}


template<class Types>
void ParserCKYAll_Impl<Types>::compute_inside_outside_probabilities()
{
  compute_inside_probabilities();
  static int start_symbol = SymbolTable::instance_nt().get(LorgConstants::tree_root_name);
  chart->get_root().get_edge(start_symbol).uedge().get_annotations().reset_outside_probabilities(1.0);
  compute_outside_probabilities();
}

#endif