// -*- mode: c++ -*-
#ifndef _MAXRULEPROBABILITYKB_CPP_
#define _MAXRULEPROBABILITYKB_CPP_

#include "MaxRuleProbabilityKB.h"


double MaxRuleProbabilityKB::log_normalisation_factor = 0;
unsigned MaxRuleProbabilityKB::size = 0;


std::ostream & MaxRuleProbabilityKB::operator>> (std::ostream & out) const
{
  for(auto& cand: candidates) { out << "cand:" << cand.probability << " "; }
  return out;
}


void MaxRuleProbabilityKB::update_lexical(LBEdge& e, const LexicalDaughter& dtr)
{
   //BLOCKTIMING("MaxRuleProbabilityKB::update_lexical");
 const AnnotationInfo & a = e.get_annotations();
  const LexicalRuleC2f* rule = dtr.get_rule();
  assert(rule != NULL);

  packed_edge_probability_with_index pep;
  pep.probability = QInsideComputer::compute(a, rule, log_normalisation_factor);

//    std::cout << "lexical " << pep.probability << std::endl;
  assert(pep.probability <=0);

  pep.dtrs = &dtr;

  candidates.push_back(pep);

  if(derivations.empty())
     derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
      derivations[0] = pep;

//   std::cout << *this << std::endl;
}

void MaxRuleProbabilityKB::update_unary(UEdge& e, const UnaryDaughter& dtr)
{
  //BLOCKTIMING("MaxRuleProbabilityKB::update_unary");
  const AnnotationInfo & a = e.get_annotations();
  packed_edge_probability_with_index pep;
  pep.dtrs = &dtr;
  //  std::cout << "before ump" << std::endl;
  pep.probability= QInsideComputer::compute(a, dtr, log_normalisation_factor);

//    std::cout << "unary "<< pep.probability << std::endl;
  assert(pep.probability <=0);

  candidates.push_back(pep);

  if(derivations.empty())
       derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
       derivations[0] = pep;
  
//   std::cout << *this << std::endl;
}

void MaxRuleProbabilityKB::update_binary(LBEdge& e, const BinaryDaughter& dtr)
{
  //BLOCKTIMING("MaxRuleProbabilityKB::update_binary");
  const AnnotationInfo & a = e.get_annotations();
  packed_edge_probability_with_index pep;
  pep.dtrs = &dtr;

  pep.probability= QInsideComputer::compute(a, dtr, log_normalisation_factor);

  //  std::cout << candidates.size() << std::endl;
//   std::cout << "binary " << pep.probability << std::endl;
  //  std::cout << pep.dtrs << std::endl;


  //  std::cout << pep.probability << std::endl;
  assert(pep.probability <=0);

  candidates.push_back(pep);


  if(derivations.empty())
    derivations.push_back(pep);
  else if(pep.probability > derivations[0].probability)
    derivations[0] = pep;

//     std::cout << *this << std::endl;
}

struct gt_pep
{
  bool operator()(const packed_edge_probability_with_index& p1, const packed_edge_probability_with_index& p2) const
  {
    return p1 > p2;
  }
};



void MaxRuleProbabilityKB:: finalize()
{
  //  std::cout << "size candidates: " << candidates.size() << std::endl;

  if(!candidates.empty()) {
    if(candidates.size() > size) {

      // std::cout << "BEFORE NTH " << size << std::endl;
      // for (unsigned i = 0; i < candidates.size(); ++i)
      //   {
      //     std::cout << candidates[i].probability << " ";

      //     if(candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(candidates[i].dtrs->is_binary())
      //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(!candidates[i].dtrs->is_binary() && !candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     std::cout << std::endl;
      //   }

      std::nth_element(candidates.begin(),candidates.begin()+size,candidates.end(), gt_pep());
      candidates.resize(size);

      //       std::cout << "AFTER NTH " << size << std::endl;
      // for (int i = 0; i < candidates.size(); ++i)
      //   {
      //     std::cout << candidates[i].probability << " ";

      //     if(candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const LexicalPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(candidates[i].dtrs->is_binary())
      //       std::cout << *(static_cast<const BinaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     if(!candidates[i].dtrs->is_binary() && !candidates[i].dtrs->is_lexical())
      //       std::cout << *(static_cast<const UnaryPackedEdgeDaughters*>(candidates[i].dtrs)->get_rule());
      //     std::cout << std::endl;
      //   }



    }
    std::make_heap(candidates.begin(),candidates.end());

    std::pop_heap(candidates.begin(),candidates.end());
    candidates.pop_back();
  }
}

void MaxRuleProbabilityKB::extend_derivation(PEdge* edge, unsigned i, bool licence_unaries)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability_with_index& last = derivations[derivations.size() -1];

    //    std::cout << "last.probability " << last.probability << std::endl;

    assert(last.probability <= 0);

    find_succ(edge,last,licence_unaries);
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

void MaxRuleProbabilityKB::find_succ(PEdge* edge, packed_edge_probability_with_index& pep, bool licence_unaries)
{
  if(pep.dtrs->is_lexical())  { return;}
  // binary -> extend left and right daughters

  if(pep.dtrs->is_binary()) {
    const BinaryDaughter* d = static_cast<const BinaryDaughter*>(pep.dtrs);

    //extend to the left
    PEdge& left  = d->left_daughter();
    unsigned nextleft = pep.left_index + 1;
    left.extend_derivation(nextleft+1,true);

    // we haven't reached the expected number of solutions
    if(nextleft < left.get_best().n_deriv()) {

      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = QInsideComputer::compute(static_cast<LBEdge *>(edge)->get_annotations(), *d, log_normalisation_factor, p.left_index, p.right_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }

    }

    //extend to the right
    PEdge& right = d->right_daughter();
    unsigned nextright = pep.right_index + 1;

    right.extend_derivation(nextright+1,true);

    if(nextright < right.get_best().n_deriv()) {
      //        std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability_with_index p(pep);
      p.right_index = nextright;
      p.probability = QInsideComputer::compute(static_cast<LBEdge *>(edge)->get_annotations(), *d, log_normalisation_factor, p.left_index, p.right_index);

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
    LBEdge& left  = d->left_daughter();
    unsigned nextleft = pep.left_index + 1;

    left.extend_derivation(nextleft+1, false);

    if(nextleft < left.get_best().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = QInsideComputer::compute(static_cast<UEdge *>(edge)->get_annotations(), *d, log_normalisation_factor, p.left_index);

      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;


      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }
}



#endif /* _MAXRULEPROBABILITYKB_CPP_ */
