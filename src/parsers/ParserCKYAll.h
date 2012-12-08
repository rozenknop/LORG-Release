// -*- mode: c++ -*-
#ifndef PARSERCKYALL_H
#define PARSERCKYALL_H

#include "grammars/GrammarAnnotated.h"
#include "ParserCKY.h"


#include "ChartCKY.h"

#include "edges/PackedEdge.h"


#include "rules/BRule.h"
#include "rules/URule.h"

#include "utils/RuleVect2Map.h"

#include "utils/data_parsers/BURuleInputParser.h"

#include <algorithm>
#include <numeric>

//#define USE_THREADS 1

#ifdef USE_THREADS
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
using namespace tbb;
#endif


typedef std::vector< std::vector<std::vector< std::vector<unsigned> > > > annot_descendants_type;



class ParserCKYAll : public ParserCKY< GrammarAnnotated<BRuleC2f,URuleC2f, LexicalRuleC2f> >
{
 public:
  typedef GrammarAnnotated<BRuleC2f,URuleC2f, LexicalRuleC2f> AGrammar;
  typedef ParserCKY<AGrammar> Parser;

  /**
     \brief ParserCKYAll destructor
     */
  virtual ~ParserCKYAll() {};

  /**
     \brief ParserCKYAll constructor
     \param cgs the grammars from the coarsest to the finest
     \param prior_map, the prior values for each NT
     \param beam_threshold, the threshold for prior-based pruning
     \param annot_descendants the data structures for C2F mappings
     \param accurate, triggers accurate mode
     \param min_beam, minimum length for sentences to apply pruning
     \param stubborn, number of retries
  */
  ParserCKYAll(std::vector<AGrammar*>& cgs, const std::vector<double>& prior_map, double beam_threshold,
               const annot_descendants_type& annot_descendants,
               bool accurate, unsigned min_beam, int stubborn);

  /**
     \brief parses the sentence using the grammar
     \param start_symbol the axiom symbol of the grammar
  */
  virtual void parse(int start_symbol) const = 0;


  /**
     \brief Performs the PCFG-LA coarse-to-fine beam
  */
  virtual void beam_c2f(int start_symbol) = 0;


  /**
     \brief computes the best solution (depends on the algorithm)
  */
  virtual void extract_solution()=0;

  virtual   void initialise_chart(const std::vector< Word >& s,
                                  const std::vector<bracketing>& bs) = 0;

  virtual bool is_chart_valid(int start_symbol) = 0;

  virtual void get_parses(int start_symbol, unsigned kbest, bool always_output_forms, bool output_annotations,
                          std::vector<std::pair<PtbPsTree *,double> >& best_trees) = 0;

  virtual void clean() = 0;


 protected: // attributes
  std::vector<AGrammar *> grammars; ///< the grammars to beam
  const std::vector<double> priors;  ///< used for prior-based beam

  double prior_beam_threshold; ///< the default threshold value for prior-based pruning


  annot_descendants_type annot_descendants; // fast access to annotation histories

  std::vector<double> io_beam_thresholds; ///< thresholds used for c2f IO beam (one per c2f level)

  bool accurate;  ///< triggers accurate/slow parsing thresholds

  unsigned min_length_beam; ///< minimum sentence length to apply beam on

  int stubbornness; ///< number of tries with incremental prior-based pruning
};


/**
 \class ParserCKYAll
 \brief represents a parsing device for probabilistic cfgs using the cky algorithm
 */
template<class Types>
class ParserCKYAll_Impl : public ParserCKYAll
{
public:
  typedef typename Types::BRule BinaryRule;
  typedef typename Types::URule UnaryRule;
  typedef typename Types::LRule LexicalRule;
  typedef typename Types::Edge Edge;
  typedef typename Types::Cell Cell;
  typedef typename Types::EdgeProbability ProbaModel;
  typedef typename Types::Chart Chart;

  /**
     \brief ParserCKYAll_Impl destructor
  */
  virtual ~ParserCKYAll_Impl() ;

  /**
     \brief ParserCKYAll_Impl constructor
     \param cgs the grammars from the coarsest to the finest
     \param prior_map, the prior values for each NT
     \param beam_threshold, the threshold for prior-based pruning
     \param annot_descendants the data structures for C2F mappings
     \param accurate, triggers accurate mode
     \param min_beam, minimum length for sentences to apply pruning
     \param stubborn, number of retries
  */
  ParserCKYAll_Impl(std::vector<AGrammar*>& cgs, const std::vector<double>& prior_map, double beam_threshold,
                    const annot_descendants_type& annot_descendants,
                    bool accurate, unsigned min_beam, int stubborn);

  /**
     \brief parses the sentence using the grammar
     \param start_symbol the axiom symbol of the grammar
  */
  void parse(int start_symbol) const;

  /**
     \brief computes the best solution (depends on the algorithm)
  */
  virtual void extract_solution() =0;

  /**
     \brief gets the idx^th grammar of the c2f-parser
     \param idx the rank of the grammar
     \return the idx^th grammar in the c2f-chain
  */
  AGrammar& get_grammar(unsigned idx);


  /**
     \brief gets the idx^th grammar of the c2f-parser
     \param idx the rank of the grammar
     \return the idx^th grammar in the c2f-chain
  */
  const AGrammar& get_grammar(unsigned idx) const;

  /**
     \brief Performs the PCFG-LA coarse-to-fine beam
  */
  void beam_c2f(int start_symbol);


  /**
     \brief returns the sentence score (depends on the parsing algorithm)
     \return a score
  */
  double get_sentence_probability() const;


  void initialise_chart(const std::vector< Word >& sentence, const std::vector<bracketing>& brackets)
  {
    chart = new Chart(sentence, get_nonterm_count(), brackets);
//     std::cout << *chart << std::endl;

  }

  void clean() { delete chart; chart = NULL;}


 private:
  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from binary rules)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary_internal(Cell& cell, bool isroot) const;

  /** \brief Add unary rules at this position in the chart
      (only consider non-terminals created from pos)
      \param cell the cell to fill
      \param isroot true if cell is root
  */
  void add_unary_init(Cell& cell, bool isroot) const;

  /**
     \brief processes the internal rules, adding edges to the chart where appropriate
     \param beam_threshold the prior beam threshold to prune the base packed forest
  */
  void process_internal_rules(double beam_threshold) const;
  void process_cell(Cell& cell, double beam_threshold) const;

  /**
     \brief fill the result cell with the most probable edges for each lhs,
     created from edges contained in left_cell and right_cell
     \param left_cell the  leftmost cell to combine
     \param right_cell the rightmost cell to combine
     \param result_cell the cell to store new edges
  */
  void get_candidates(Cell& left_cell,
                      Cell& right_cell,
                      Cell& result_cell) const;

  /**
     \brief process all unary rules for this edge - unary chains are not followed
     \param cell the cell to fill
     \param lhs the entry in the cell
     \param isroot true if cell is the rootcell
  */
  void process_unary(Cell& cell, int lhs, bool isroot) const;


protected:
  /**
     \brief computes the inside probability for all nodes in chart
  */
  void compute_inside_probabilities();

  /**
     \brief computes the outside probability for all nodes in chart
  */
  virtual void compute_outside_probabilities();


  /**
     \brief computes the inside & outside scores in the chart
  */
  void compute_inside_outside_probabilities();


  /**
     \brief prunes the chart according to I/O scores of nodes
     if a node I/0 is lower by k% than the best node in the cell, remove it
  */
  void beam_chart_io_relative() const;


  /**
     \brief beam the chart
     \param sent_prob, the score of the current sentence (with the base grammar)
     \param threshold, threshold for pruning
     \param huang, triggers the use of huang/charniak pruning with the base grammar
  */
  void beam_chart(double sent_prob, double threshold, bool huang);


  /**
     \brief creates a c2f grammar from a vector of grammar
     \param cgs, grammars from coarsest to finest
  */
  void create_coarse_to_fine_mapping(std::vector<AGrammar*>& cgs);


  /**
     \brief Performs the PCFG-LA coarse-to-fine beaming
     \para current_grammars the grammars from coarsest to finest
     \para current_annot_descendants the histories of annotations
  */
  void beam_c2f(const std::vector<AGrammar*>& current_grammars,
                const annot_descendants_type& current_annot_descendants);

  void change_rules_resize(unsigned step,
                           const std::vector<AGrammar*>& current_grammars) const;


  bool is_chart_valid(int start_symbol) {return chart->is_valid(start_symbol);}


  void get_parses(int start_symbol, unsigned kbest, bool always_output_forms, bool output_annotations,
                  std::vector<std::pair<PtbPsTree *,double> >& best_trees);


 protected: // attributes
  Chart * chart; // the chart
};

#endif /*PARSERCKYALL_H*/
