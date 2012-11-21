// -*- mode: c++ -*-
#ifndef PACKEDEDGE_H_
#define PACKEDEDGE_H_

#include "PackedEdgeDaughters.h"

#include "grammars/AnnotatedLabelsInfo.h"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "AnnotationInfo.h"

#include "rules/BRuleC2f.h"
#include "rules/URuleC2f.h"
#include "rules/LexicalRuleC2f.h"

#include "utils/PtbPsTree.h"

#include "PCKYAllCell.h"



#include "utils/SymbolTable.h"
#include "PackedEdgeProbability.h"

#include "utils/lorg_functional.h"

typedef std::pair< int, unsigned> asymb;
typedef boost::unordered_map<asymb, boost::unordered_map< asymb, asymb> > PathMatrix;


/**
  \class PackedEdge
  \brief represents an edge in a chart
*/



template<class PEP>
class PackedEdge
{
public:
  typedef PEP ProbaModel;
  typedef PCKYAllCell<PackedEdge<PEP> > Cell;
  typedef PackedEdge<PEP> Edge;
  typedef PackedEdgeDaughters Daughters;
  typedef BinaryPackedEdgeDaughters<Cell> BinaryDaughters;
  typedef UnaryPackedEdgeDaughters<Cell>  UnaryDaughters;
  typedef LexicalPackedEdgeDaughters LexicalDaughters;

protected:
  typedef std::vector<BinaryDaughters> bvector;
  typedef std::vector<UnaryDaughters>  uvector;
  typedef std::vector<LexicalDaughters>  lvector;

  typedef typename bvector::iterator biterator;
  typedef typename uvector::iterator  uiterator;
  typedef typename lvector::iterator  literator;

  typedef typename bvector::const_iterator const_biterator;
  typedef typename uvector::const_iterator  const_uiterator;
  typedef typename lvector::const_iterator  const_literator;

public:

  /**
     \brief Constructor for creating an empty edge
  */
  PackedEdge()
  {
    this->local_resize_annotations(1);
  }

  /**
     \brief Constructor
     \param ped daughters to add
  */
    PackedEdge(const BinaryDaughters& ped)
  {
    binary_daughters.push_back(ped);
    this->local_resize_annotations(1);
  }

  /**
     \brief Constructor
     \param ped daughters to add
  */
  PackedEdge(const UnaryDaughters& ped)
  {
    unary_daughters.push_back(ped);
    this->local_resize_annotations(1);
  }

  /**
     \brief Constructor
     \param ped daughters to add
  */
  PackedEdge(const LexicalDaughters& ped)
  {
    lexical_daughters.push_back(ped);
    this->local_resize_annotations(1);
  }

  /**
     \brief Destructor
  */
  ~PackedEdge()  {}


  /**
     \brief get the daughters of the edge
  */
  const bvector& get_binary_daughters() const;
  bvector& get_binary_daughters();

  const uvector& get_unary_daughters() const;
  uvector& get_unary_daughters();

  const lvector& get_lexical_daughters() const;
  lvector& get_lexical_daughters();



 /**
     \brief get one particular daughter of the edge
  */
  BinaryDaughters& get_binary_daughter(unsigned i);
  const BinaryDaughters& get_binary_daughter(unsigned i) const;

  UnaryDaughters& get_unary_daughter(unsigned i);
  const UnaryDaughters& get_unary_daughter(unsigned i) const;

  LexicalDaughters& get_lexical_daughter(unsigned i);
  const LexicalDaughters& get_lexical_daughter(unsigned i) const;

  /**
     \brief is the daughter lexical ?
  */
  bool get_lex() const;

  /**
     \brief resize the vector of annotations
     \param size the new size
   */
  void local_resize_annotations(unsigned size);

  /**
     \brief Output operator
     \param out the ostream to write on
     \param edge the edge object to write
     \return the used ostream
  */
  template<class OPEP>
  friend std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge);



  /**
     \brief return an AnnotationInfo (inside/outside probs +scaling values)
   */
  AnnotationInfo& get_annotations();
  const AnnotationInfo& get_annotations() const;

  /**
     \brief get previously calculated annotations (maxn parsing)
   */
  std::vector<AnnotationInfo>& get_annotations_backup();
  const std::vector<AnnotationInfo>& get_annotations_backup() const;


  /**
     \brief store annotions for later use (maxn parsing)
   */
  void backup_annotations()
  {
    best.backup_annotations(annotations);
  };


  /**
     \brief reset annotation probabilities at the current node
     \param value the value to reset annotations to
  */
  void reset_probabilities(const double& value);

  /**
     \brief reset annotation probabilities at the current node
     \param value the value to reset annotations to
   */
  void reset_inside_probabilities(const double& value);

  /**
     \brief reset annotation probabilities at the current node
     \param value the value to reset annotations to
  */
  void reset_outside_probabilities(const double& value);





  /**
     \brief computes the inside probability of a packed forest
  */
  void prepare_inside_probability();
  void adjust_inside_probability();


  /**
     \brief computes the outside probability of a packed forest
  */
  void prepare_outside_probability();
  void adjust_outside_probability();

  /**
     \brief build and add a daughter (binary, unary and lexical versions)
   */
  void add_daughters(Cell * left,
                     Cell * right, const BRuleC2f* rule)
  {
    binary_daughters.push_back(BinaryDaughters(left,right,rule));
  }

  void add_daughters(Cell * left, const URuleC2f* rule)
  {
    unary_daughters.push_back(UnaryDaughters(left,rule));
  }

  void add_daughters(const LexicalRuleC2f* rule, const Word* w)
  {
    lexical_daughters.push_back(LexicalDaughters(rule, w));
  }

  /**
     \brief get the structure holding the "best calculation for the prob. model" whatever it means
   */
  const PEP& get_prob_model() const;
  PEP& get_prob_model();


  /**
     \brief set unary chains
     \note TODO remove "viterbi" from name
   */
  static void set_viterbi_unary_chains(const PathMatrix& pathmatrix);
  static const PathMatrix& get_viterbi_unary_chains();


  void replace_rule_probabilities(unsigned i);

  void create_viterbi(unsigned size);
  void create_maxrule();
  void create_maxrule_kb();
  void create_maxrule_mult();

  void find_best_multiple_grammars();
  void find_best_multiple_grammars_lexical();
  void find_best_multiple_grammars_binary();
  void find_best_multiple_grammars_unary();

  void compute_best_lexical();
  void compute_best_unary();
  void compute_best_binary();


  void extend_derivation(unsigned i, bool licence_unaries)
  {
    best.extend_derivation(this,i, licence_unaries);
  }

  bool valid_prob_at(unsigned i) const;

  void clean_invalidated_binaries();
  
  PtbPsTree * to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const;

  bool has_solution(unsigned i) const ;


protected :
  bvector binary_daughters;    ///< set of possible daughters

  uvector unary_daughters;     ///< set of possible daughters
  lvector lexical_daughters;   ///< a possible lexical daughter
  AnnotationInfo annotations;  ///< probabilities

  static PathMatrix viterbi_unary_chains;


private:
  PEP best;

  void to_ptbpstree(PtbPsTree& tree, PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
		    bool append_annot, bool outpu_forms) const;

public:
  void process(function<void(Edge &, BinaryDaughters &)> f) { for(auto& d: get_binary_daughters()) f(*this, d); }
  void process(function<void(Edge &, UnaryDaughters &)> f) { for(auto& d: get_unary_daughters()) f(*this, d); }
  void process(function<void(Edge &, LexicalDaughters &)> f) {for(auto& d: get_lexical_daughters()) f(*this, d);}

  void process(function<void(BinaryDaughters &, AnnotationInfo &)> f) { for(auto& d: get_binary_daughters()) f(d, get_annotations()); }
  void process(function<void(UnaryDaughters &, AnnotationInfo &)> f) { for(auto& d: get_unary_daughters()) f(d, get_annotations()); }
  void process(function<void(LexicalDaughters &, AnnotationInfo &)> f) {for(auto& d: get_lexical_daughters()) f(d, get_annotations());}
  
  void process(function<void(PEP &, Edge &, const LexicalDaughters &)> f) {for(auto& d: get_lexical_daughters()) f(get_prob_model(), *this, d);}
  void process(function<void(PEP &, Edge &, const BinaryDaughters &)> f) {for(auto& d: get_binary_daughters()) f(get_prob_model(), *this, d);}
  void process(function<void(PEP &, Edge &, const UnaryDaughters &)> f) {for(auto& d: get_unary_daughters()) f(get_prob_model(), *this, d);}
  void process(function<void(PEP &)> f) {f(get_prob_model());}

  void process(function<void(Edge &)> f) { f(*this); }

  template<typename Function, typename... OtherFunctions>
  void apply(Function&& f, OtherFunctions&&... o) {process(toFunc(f));apply(o...);}
  void apply() {}
};




template<class PEP>
inline
AnnotationInfo& PackedEdge<PEP>::get_annotations() {return annotations;}

template<class PEP>
inline
const AnnotationInfo& PackedEdge<PEP>::get_annotations() const {return annotations;}

template<class PEP>
inline
std::vector<AnnotationInfo>& PackedEdge<PEP>::get_annotations_backup() {return best.get_annotations_backup();}

template<class PEP>
inline
const std::vector<AnnotationInfo>& PackedEdge<PEP>::get_annotations_backup() const {return best.get_annotations_backup();}

template<class PEP>
inline
const typename PackedEdge<PEP>::bvector& PackedEdge<PEP>::get_binary_daughters() const {return binary_daughters;}


template<class PEP>
inline
typename PackedEdge<PEP>::bvector& PackedEdge<PEP>::get_binary_daughters() {return binary_daughters;}

template<class PEP>
inline
const typename PackedEdge<PEP>::uvector& PackedEdge<PEP>::get_unary_daughters() const {return unary_daughters;}

template<class PEP>
inline
typename PackedEdge<PEP>::uvector& PackedEdge<PEP>::get_unary_daughters() {return unary_daughters;}

template<class PEP>
inline
typename PackedEdge<PEP>::BinaryDaughters& PackedEdge<PEP>::get_binary_daughter(unsigned i) {return binary_daughters[i];}

template<class PEP>
inline
const typename PackedEdge<PEP>::BinaryDaughters& PackedEdge<PEP>::get_binary_daughter(unsigned i) const {return binary_daughters[i];}


template<class PEP>
inline
typename PackedEdge<PEP>::UnaryDaughters& PackedEdge<PEP>::get_unary_daughter(unsigned i) {return unary_daughters[i];}

template<class PEP>
inline
const typename PackedEdge<PEP>::UnaryDaughters& PackedEdge<PEP>::get_unary_daughter(unsigned i) const {return unary_daughters[i];}

template<class PEP>
inline
bool PackedEdge<PEP>::get_lex() const {return !(lexical_daughters.empty());}

template<class PEP>
inline
const typename PackedEdge<PEP>::lvector& PackedEdge<PEP>::get_lexical_daughters() const {return lexical_daughters;}

template<class PEP>
inline
typename PackedEdge<PEP>::lvector& PackedEdge<PEP>::get_lexical_daughters() {return lexical_daughters;}

template<class PEP>
inline
typename PackedEdge<PEP>::LexicalDaughters& PackedEdge<PEP>::get_lexical_daughter(unsigned i) {return lexical_daughters[i];}

template<class PEP>
inline
const typename PackedEdge<PEP>::LexicalDaughters& PackedEdge<PEP>::get_lexical_daughter(unsigned i) const {return lexical_daughters[i];}

template<class PEP>
inline
void PackedEdge<PEP>::local_resize_annotations(unsigned size) {annotations.resize(size);}

/////////////////////////////////////////////////////
// parsing
////////////////////////////////////////////////

template<class PEP>
inline
const PEP& PackedEdge<PEP>::get_prob_model() const {return best;}

template<class PEP>
inline
PEP& PackedEdge<PEP>::get_prob_model() {return best;}

template<class PEP>
inline
bool PackedEdge<PEP>::valid_prob_at(unsigned i) const
{
  return get_annotations().valid_prob_at(i, LorgConstants::NullProba);
}


//////////

template <class PEP>
PathMatrix PackedEdge<PEP>::viterbi_unary_chains = PathMatrix();

template <class PEP>
void PackedEdge<PEP>::reset_inside_probabilities(const double& value)
{
  annotations.reset_inside_probabilities(value);
}

template <class PEP>
void PackedEdge<PEP>::reset_outside_probabilities(const double& value)
{
  annotations.reset_outside_probabilities(value);
}

template <class PEP>
void PackedEdge<PEP>::reset_probabilities(const double& value)
{
  annotations.reset_inside_probabilities(value);
  annotations.reset_outside_probabilities(value);
}




template <class PEP>
void PackedEdge<PEP>::prepare_inside_probability()
{
  this->get_annotations().inside_probabilities_unary_temp.array = this->get_annotations().inside_probabilities.array ;
  for(auto & prob: this->get_annotations().inside_probabilities_unary_temp.array) {
    if (prob != LorgConstants::NullProba) prob = 0;
  }
}

template <class PEP>
void PackedEdge<PEP>::adjust_inside_probability()
{
  for (unsigned i = 0; i < this->get_annotations().inside_probabilities.array.size(); ++i)
    {
      if(this->get_annotations().inside_probabilities.array[i] != LorgConstants::NullProba)
        this->get_annotations().inside_probabilities.array[i] += this->get_annotations().inside_probabilities_unary_temp.array[i];
    }
}


template <class PEP>
void PackedEdge<PEP>::prepare_outside_probability()
{
  this->get_annotations().outside_probabilities_unary_temp.array = this->get_annotations().outside_probabilities.array ;
  for(auto & prob: this->get_annotations().outside_probabilities_unary_temp.array) {
    if (prob != LorgConstants::NullProba) prob = 0;
  }
}

template <class PEP>
void PackedEdge<PEP>::adjust_outside_probability()
{
  for (unsigned i = 0; i < this->get_annotations().outside_probabilities.array.size(); ++i)
    {
      if(this->get_annotations().outside_probabilities.array[i] != LorgConstants::NullProba)
        this->get_annotations().outside_probabilities.array[i] += this->get_annotations().outside_probabilities_unary_temp.array[i];
    }
}





// should be renamed (Why Viterbi?) and moved
template <class PEP>
void PackedEdge<PEP>::set_viterbi_unary_chains(const PathMatrix& pathmatrix)
{
  viterbi_unary_chains = pathmatrix;
}

template <class PEP>
const PathMatrix& PackedEdge<PEP>::get_viterbi_unary_chains()
{
  return viterbi_unary_chains;
}


struct c2f_replace_struct_helper
{
  unsigned idx;
  c2f_replace_struct_helper(unsigned i) : idx(i) {};

  //reset the rule to the 'finer_id'th rule in the vector
  //TODO implement the Charniak-style coarse-to-fine method - this will involve replacing this rule with several
  template <typename T>
  void operator ()(T& daughter) const
  {
    daughter.set_rule(daughter.get_rule()->get(idx));
  }
};

template <class PEP>
void PackedEdge<PEP>::replace_rule_probabilities(unsigned i)
{
  // for all possible daughters
  c2f_replace_struct_helper c2f_replacer(i);
  std::for_each(binary_daughters.begin(),binary_daughters.end(), c2f_replacer);
  std::for_each(unary_daughters.begin(),unary_daughters.end(),   c2f_replacer);
  std::for_each(lexical_daughters.begin(), lexical_daughters.end(),   c2f_replacer);
}

template <class PEP>
void PackedEdge<PEP>::create_viterbi(unsigned size)
{
  best.set_size(size);
}

template <class PEP>
void PackedEdge<PEP>::create_maxrule()
{}

template <class PEP>
void PackedEdge<PEP>::create_maxrule_kb()
{}

template <class PEP>
void PackedEdge<PEP>::create_maxrule_mult()
{}

template <class PEP>
void PackedEdge<PEP>::find_best_multiple_grammars()
{
  best.pick_best(this);
}

template <class PEP>
void PackedEdge<PEP>::find_best_multiple_grammars_lexical()
{
  for(literator it(lexical_daughters.begin()); it != lexical_daughters.end(); ++it) {
    best.pick_best_lexical(this, *it);
  }
}

template <class PEP>
void PackedEdge<PEP>::find_best_multiple_grammars_binary()
{
  for(biterator it(binary_daughters.begin()); it != binary_daughters.end(); ++it) {
    // calculate the probability for each edge in this packed edge
    // if it's the best probability so far, update the best edge info
    best.pick_best_binary(this,*it);
  }
}

template <class PEP>
void PackedEdge<PEP>::find_best_multiple_grammars_unary()
{
  //  std::cout << "before updates" << std::endl;
  //  std::cout << "best: " << best << std::endl;

  for(uiterator it(unary_daughters.begin()); it != unary_daughters.end(); ++it) {
    best.pick_best_unary(this,*it);
  }

  best.pick_best();
}




template <class PEP>
void PackedEdge<PEP>::compute_best_lexical()
{
  for(literator it(lexical_daughters.begin()); it != lexical_daughters.end(); ++it) {
    best.update(this->get_annotations(), *it);
  }
}

template <class PEP>
void PackedEdge<PEP>::compute_best_unary()
{
  //  std::cout << "before updates" << std::endl;
  //  std::cout << "best: " << best << std::endl;

  for(uiterator it(unary_daughters.begin()); it != unary_daughters.end(); ++it) {
    best.update(this->get_annotations(),*it);
  }
  //  std::cout << "before finalize" << std::endl;
  //  if(best)
  //best.finalize(this);
  best.finalize();
}

template <class PEP>
void PackedEdge<PEP>::compute_best_binary()
{
  for(biterator it(binary_daughters.begin()); it != binary_daughters.end(); ++it) {
    // calculate the probability for each edge in this packed edge
    // if it's the best probability so far, update the best edge info
    best.update(this->get_annotations(),*it);
  }
}

/////////////
/// clean
////////////

template <class PEP>
void PackedEdge<PEP>::clean_invalidated_binaries()
{
  binary_daughters.erase(std::remove_if(binary_daughters.begin(), binary_daughters.end(), toFunc(& BinaryDaughters::points_towards_invalid_cells)), binary_daughters.end());
  
  
  // Reclaim memory !
  if(binary_daughters.capacity() != binary_daughters.size()) {
    decltype(binary_daughters) tmp;
    tmp.swap(binary_daughters);
    binary_daughters.insert(binary_daughters.begin(), tmp.begin(), tmp.end());
  }
}




///////////////////////////////
////// printing
//////////////////////////////

#include <sstream>


// decode the path-matrix to read a unary chain between start and stop
// TODO: create a class for the matrix and encapsulate
// (because matrix for max-rule has a different type)
inline
unsigned decode_path(PtbPsTree& tree,
		     PtbPsTree::depth_first_iterator& pos,
		     const PathMatrix& paths,
		     const asymb& start,
		     const asymb& end,
		     bool append_annot
		     )
{
  bool final = false;
  asymb candidate = start;
  unsigned path_length = 0;


  while(!final) {

    boost::unordered_map<asymb, boost::unordered_map< asymb, asymb> >::const_iterator it1 =
      paths.find(candidate);

    if(it1 == paths.end())
      final = true;

    else {
      boost::unordered_map< asymb, asymb>::const_iterator it2 = it1->second.find(end);

      if(it2 == it1->second.end())
	final = true;
      else {
	//	std::cout << "adding node " << SymbolTable::instance_nt()->translate(it2->second.first) << std::endl;

	std::ostringstream node_content;
	node_content << SymbolTable::instance_nt().translate(it2->second.first);
	if(append_annot) node_content << "_" << it2->second.second;

	pos=tree.add_last_daughter(pos,node_content.str());
	candidate = it2->second;
	++path_length;
      }
    }
  }

  return path_length;

}


// to prevent a bug on maia
// it seems that inf and nan are the same on the cluster :(
bool my_isinvalid(const double& input)
{

    bool b1 = std::isnan(input);
    bool b2 = input != input;

    bool b3 = false;
    if(- std::numeric_limits<double>::infinity() == input) {
      //    std::cout << "-inf" << std::endl;
        b3 = true;
    }

    if(std::numeric_limits<double>::infinity() == input) {
      //        std::cout << "+inf" << std::endl;
	b3 = false;
    }

    return b3 || b1 || b2;


}


template <class PEP>
PtbPsTree * PackedEdge<PEP>::to_ptbpstree(int lhs, unsigned ith_deriv, bool append_annot, bool output_forms) const
{

  PtbPsTree * tree = NULL;

  std::ostringstream node_content;
  node_content << SymbolTable::instance_nt().translate(lhs) ;
  // we assume that root is TOP hence one annotation: 0
  if(append_annot) node_content << '_' << 0;

  //  std::cout << SymbolTable::instance_nt()->translate(get_lhs()) << std::endl;
  //  std::cout << "size daughters: " << best_dtrs_vector.size() << std::endl;
  //  if(best_dtrs_vector[0] == NULL) std::cout << "daughter is NULL" << std::endl;


  //  std::cout << "log prob deriv: " << best.get(ith_deriv).probability << std::endl;

  if(best.get(ith_deriv).probability == - std::numeric_limits<double>::infinity())
    return NULL;

  if(best.get(ith_deriv).dtrs == NULL) {
    //    std::cerr << "no daughters" << std::endl;
    return NULL;
  }
  if(std::isnan(best.get(ith_deriv).probability)) {
    //    std::cerr << "invalid prob" << std::endl ;
    return NULL;
  }
  // if(my_isinvalid(best.get(ith_deriv).probability)) {
  //   //g    std::cerr << "invalid prob" << std::endl;
  //   return NULL;
  // }


  //  std::cout << best.get(ith_deriv).probability << std::endl;


    tree = new PtbPsTree(node_content.str());
    PtbPsTree::depth_first_iterator pos = tree->dfbegin();

    //    std::cout << "allocation done" << std::endl;

    if(best.get(ith_deriv).dtrs->is_binary()) {

      // std::cout << "deriv " << ith_deriv
      //           << "\t" << best.get(ith_deriv).get_left_index()
      //           << "\t" << best.get(ith_deriv).get_right_index() << std::endl;

      const BinaryDaughters * daughters =  static_cast<const BinaryDaughters*>(best.get(ith_deriv).dtrs);
      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter()->get_edge(rhs0).to_ptbpstree(*tree, pos, rhs0, best.get(ith_deriv).get_left_index(), append_annot, output_forms);
      int rhs1 = daughters->get_rule()->get_rhs1();
      daughters->right_daughter()->get_edge(rhs1).to_ptbpstree(*tree, pos, rhs1, best.get(ith_deriv).get_right_index(), append_annot, output_forms);
    }
    else {

            // std::cout << "deriv " << ith_deriv
            //     << "\t" << best.get(ith_deriv).get_left_index() << std::endl;

      const UnaryDaughters * daughters =  static_cast<const UnaryDaughters*>(best.get(ith_deriv).dtrs);
      decode_path(*tree,pos,
                  get_viterbi_unary_chains(),
                  std::make_pair(lhs,0),
                  std::make_pair(daughters->get_rule()->get_rhs0(), best.get(ith_deriv).get_left_index()),
                  append_annot);
      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter()->get_edge(rhs0).to_ptbpstree(*tree, pos, rhs0, best.get(ith_deriv).get_left_index(), append_annot, output_forms);
    }

    return tree;
}

template <class PEP>
void PackedEdge<PEP>::to_ptbpstree(PtbPsTree& tree,
			      PtbPsTree::depth_first_iterator& pos, int lhs, unsigned index,
			      bool append_annot, bool output_forms) const
{
  std::ostringstream node_content;

  // default height for unary chains
  unsigned added_height = 1;

  // either a valid internal node or a lexical one
  assert(best.get(index).dtrs || get_lex());

  // if we can't find a branching (lexical node)
  // the nwe exit
  if(!best.get(index).dtrs) {return;}


  // if the *branching* leads to a lexical node
  // the we write the preterminal *and* the word
  if(best.get(index).dtrs->is_lexical()) {

    // std::cout << "lex: " << std::endl;
    // std::cout << "lex: " << get_lhs() << std::endl;
    // std::cout << SymbolTable::instance_nt()->translate(get_lhs()) << std::endl;
    // pos = tree.add_last_daughter(pos, SymbolTable::instance_word()->translate(get_lhs()));

    const LexicalDaughters * daughters =  static_cast<const LexicalDaughters*>(best.get(index).dtrs);


    node_content << SymbolTable::instance_nt().translate(lhs);
    if(append_annot) node_content << "_" << index;
    pos = tree.add_last_daughter(pos, node_content.str());

    const Word& w = *(daughters->get_word());
    std::string s = output_forms ? w.get_form() :
      (w.get_id() != -1 ? SymbolTable::instance_word().get_label_string(w.get_id()) : LorgConstants::token_unknown);

    pos = tree.add_last_daughter(pos, s);

    // we added two nodes here
    added_height = 2;

  } // if(best.get(index).dtrs->is_lexical())
  else {
    // std::cout << "int: " << lhs << std::endl;
    // std::cout << SymbolTable::instance_nt().translate(lhs) << std::endl;

    node_content << SymbolTable::instance_nt().translate(lhs);
    if(append_annot) node_content << "_" << index;

    pos = tree.add_last_daughter(pos, node_content.str());

    //assert(best);
    assert(best.get(index).dtrs);

    if(best.get(index).dtrs->is_binary()) {

      // std::cout << "binary" << std::endl;
      // std::cout << "deriv " << index
      //           << "\t" << best.get(index).get_left_index()
      //           << "\t" << best.get(index).get_right_index() << std::endl;

      const BinaryDaughters * daughters =  static_cast<const BinaryDaughters*>(best.get(index).dtrs);

      int rhs0 = daughters->get_rule()->get_rhs0();
      daughters->left_daughter()->get_edge(rhs0).to_ptbpstree(tree, pos, rhs0,
							best.get(index).get_left_index(), append_annot, output_forms);
      int rhs1 = daughters->get_rule()->get_rhs1();
      daughters->right_daughter()->get_edge(rhs1).to_ptbpstree(tree, pos, rhs1,
                                                         best.get(index).get_right_index(), append_annot, output_forms);
    } // if(best.get(index).dtrs->is_binary())
    else { //unary branching

      // std::cout << "unary" << std::endl;
      // std::cout << "deriv " << index
      //           << "\t" << best.get(index).get_left_index()
      //           << "\t" << best.get(index).get_right_index() << std::endl;

      const UnaryDaughters * daughters =  static_cast<const UnaryDaughters*>(best.get(index).dtrs);

      std::pair<int,int> fro = std::make_pair(lhs, index);
      std::pair<int,int> to =  std::make_pair(daughters->get_rule()->get_rhs0(),
					      best.get(index).get_left_index());

      added_height += decode_path(tree,pos,
      				  get_viterbi_unary_chains(),
      				  fro,
      				  to,
        			  append_annot);

      // std::cout << index << "\t" << best_left_indices.size() << std::endl;
      // std::cout << best_left_indices[index] << std::endl;

      int rhs0 = daughters->get_rule()->get_rhs0();
      //      std::cout << *daughters->get_rule() << std::endl;
      daughters->left_daughter()->get_edge(rhs0).to_ptbpstree(tree, pos, rhs0,
							best.get(index).get_left_index(), append_annot, output_forms);
    }

  }

  // go up because you processed the last daughter
  while(added_height--)
    pos.up();
}


// TODO refile this method above

template <class PEP>
inline
bool PackedEdge<PEP>::has_solution(unsigned i) const
{
  return best.has_solution(i);
}

template<class OPEP>
std::ostream& operator<<(std::ostream& out, const PackedEdge<OPEP>& edge)
{
  return out << edge.get_annotations().get_inside(0) << " (" << edge.best << ") ";
}

#endif /*PACKEDEDGE_H_*/
