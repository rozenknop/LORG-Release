// -*- mode: c++ -*-

#ifndef _MAXRULEMULTIPLEPROBABILITY_H_
#define _MAXRULEMULTIPLEPROBABILITY_H_

#include "PackedEdgeProbability.h"
#include "PackedEdge.h"
#include "MaxRuleTreeLogProbaComputer.h"
#include "emptystruct.h"
#include "ChartCKY.h"

// #include <numeric>
#include <vector>
#include <unordered_map>



class MaxRuleProbabilityMultiple;

struct MaxRuleMultipleTypes {
  typedef MaxRuleProbabilityMultiple Best ;
  typedef emptystruct EdgeProbability ;
  typedef emptystruct EdgeDaughterProbability ;
  typedef Word ChartWord ;
  
  typedef BRuleC2f BRule;
  typedef URuleC2f URule;
  typedef LexicalRuleC2f LRule;
  typedef PackedEdge< MaxRuleMultipleTypes > Edge ;
  typedef AnnotatedEdge< MaxRuleMultipleTypes > AEdge ;
  typedef BasePackedEdge< MaxRuleMultipleTypes > PEdge ;
  typedef UPackedEdge< MaxRuleMultipleTypes > UEdge ;
  typedef LBPackedEdge< MaxRuleMultipleTypes > LBEdge ;
  typedef PCKYAllCell< MaxRuleMultipleTypes > Cell ;
  typedef ChartCKY< MaxRuleMultipleTypes > Chart ;
  typedef BinaryPackedEdgeDaughters<MaxRuleMultipleTypes> BinaryDaughter;
  typedef UnaryPackedEdgeDaughters<MaxRuleMultipleTypes>  UnaryDaughter;
  typedef LexicalPackedEdgeDaughters<MaxRuleMultipleTypes> LexicalDaughter;
};



class MaxRuleProbabilityMultiple
{
private:
  typedef std::unordered_map<const PackedEdgeDaughters*,double> score_map_type;
  typedef std::unordered_map<const PackedEdgeDaughters*,unsigned> occ_map_type;
  typedef std::vector<packed_edge_probability_with_index> heap_type;
  typedef MaxRuleTreeLogProbaComputer<MaxRuleProbabilityMultiple> QInsideComputer;

  score_map_type scores;
  occ_map_type occ;

  std::vector<AnnotationInfo> annotations_backup;

  heap_type candidates;
  heap_type derivations;

  static double log_normalisation_factor;
  static std::vector<double> log_normalisation_factor_backup;
  static unsigned size;
  static unsigned nb_grammars;



public:
  typedef typename MaxRuleMultipleTypes::Edge Edge;
  typedef typename MaxRuleMultipleTypes::AEdge AEdge;
  typedef typename MaxRuleMultipleTypes::PEdge PEdge;
  typedef typename MaxRuleMultipleTypes::UEdge UEdge;
  typedef typename MaxRuleMultipleTypes::LBEdge LBEdge;
  typedef typename MaxRuleMultipleTypes::Cell Cell;
  typedef typename MaxRuleMultipleTypes::UnaryDaughter UnaryDaughter;
  typedef typename MaxRuleMultipleTypes::BinaryDaughter BinaryDaughter;
  typedef typename MaxRuleMultipleTypes::LexicalDaughter LexicalDaughter;
  
  MaxRuleProbabilityMultiple() : candidates(),
                                 derivations(heap_type(1))
  {candidates.reserve(50);};
  ~MaxRuleProbabilityMultiple() {};

  static void set_log_normalisation_factor(double lnf);
  static void reset_log_normalisation_factor();
  static const double& get_log_normalisation_factor();
  inline static const double& get_log_normalisation_factor(unsigned i);

  inline static void set_size(unsigned k)       {size = k;}

  inline static void set_nbgrammars(unsigned n) {nb_grammars = n;}

  inline const packed_edge_probability_with_index& get(unsigned idx) const {return derivations[idx];}
  inline       packed_edge_probability_with_index& get(unsigned idx)       {return derivations[idx];}

  inline void update_lexical(LBEdge& e, const LexicalDaughter& dtr);
  inline void update_unary(UEdge& e, const UnaryDaughter& dtr);
  inline void update_binary(LBEdge& e, const BinaryDaughter& dtr);
  inline void finalize();

  inline void pick_best_lexical(const LexicalDaughter& dtr);
  inline void pick_best_binary(const BinaryDaughter& dtr);
  inline void pick_best_unary(const UnaryDaughter& dtr);
  inline void pick_best();

  inline void find_succ(PEdge*,packed_edge_probability_with_index& pep, bool licence_unaries);
  inline void extend_derivation(PEdge*, unsigned, bool);


  inline unsigned n_deriv() const {return derivations.size();}
  inline bool has_solution(unsigned i) const
  {
    //    std::cout << "i " << i << std::endl;
    //    std::cout << "derivations.size() " << derivations.size() << std::endl;

    // if(i < derivations.size())
    //   std::cout << derivations[i].probability << std::endl;;

    return
      i < derivations.size() ;
    //&& derivations[i].probability != -std::numeric_limits<double>::infinity();
  }

  inline       std::vector<AnnotationInfo>& get_annotations_backup();
  inline const std::vector<AnnotationInfo>& get_annotations_backup() const;

  inline void backup_annotations(const AnnotationInfo& annotations);

private:
  inline void write_scores(const PackedEdgeDaughters& dtr, double probability);


  struct test_helper
  {
    const packed_edge_probability_with_index& pep;
    test_helper(const packed_edge_probability_with_index& p) : pep(p) {};

    bool operator()(const packed_edge_probability_with_index& p)
    {
      //      return false;
      return p.probability == pep.probability
        //  || p.dtrs == pep.dtrs
        ;}
  };

};
std::ostream& operator<<(std::ostream& out, const MaxRuleProbabilityMultiple & prob);



#endif /* _MAXRULEMULTIPLEPROBABILITY_H_ */
