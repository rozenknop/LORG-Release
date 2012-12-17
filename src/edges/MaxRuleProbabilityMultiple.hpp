// -*- mode: c++ -*-

#ifndef _MAXRULEMULTIPLEPROBABILITY_HPP_
#define _MAXRULEMULTIPLEPROBABILITY_HPP_

#include "MaxRuleProbabilityMultiple.h"
#include "edges/PackedEdge.hpp"

const double& MaxRuleProbabilityMultiple::get_log_normalisation_factor(unsigned i)
{
  //  std::cout << "size: " << log_normalisation_factor_backup.size() << std::endl;;
  return log_normalisation_factor_backup[i];
}

void MaxRuleProbabilityMultiple::write_scores(const PackedEdgeDaughters& dtr, double probability)
{
  //  if (probability == -std::numeric_limits<double>::infinity())
  //    return;

  std::pair<score_map_type::iterator, bool> res = scores.insert(std::make_pair(&dtr,probability));

  if(!res.second) { // dtr is already in the map -> add probability
    //    std::cout << "previous is " << (res.first)->second ;
    (res.first)->second += probability;
    //    std::cout << "new is " << (res.first)->second  << std::endl;
  }

  occ[&dtr]++;

}


void MaxRuleProbabilityMultiple::update_lexical(LBEdge& e, const LexicalDaughter& dtr)
{
  const AnnotationInfo & a = e.get_annotations();
  double probability(QInsideComputer::compute(a, dtr.get_rule(), log_normalisation_factor));

  if (probability > derivations[0].probability) {
    derivations[0].probability = probability;
    derivations[0].dtrs = &dtr;
  }

  write_scores(dtr, probability);
}

void MaxRuleProbabilityMultiple::update_unary(UEdge& e, const UnaryDaughter & dtr)
{
  const AnnotationInfo & a = e.get_annotations();
  double probability = -std::numeric_limits<double>::infinity();

  LBEdge& left  = dtr.lbdaughter();
  if(left.get_best().get(0).dtrs && (left.get_best().get(0).dtrs->is_lexical() || left.get_best().get(0).dtrs->is_binary())) {
    probability =  QInsideComputer::compute(a, dtr, log_normalisation_factor);
  }
  // else {
  //   std::cout << "ERROR" << std::endl;
  // }

  // if(std::isnan(probability))
  //   return;

  if (probability > derivations[0].probability) {
    derivations[0].probability = probability;
    derivations[0].dtrs = &dtr;
  }

  write_scores(dtr, probability);
}

void MaxRuleProbabilityMultiple::update_binary(LBEdge& e, const BinaryDaughter& dtr)
{
  const AnnotationInfo & a = e.get_annotations();

  //calculate the probability for this edge -
  //if it's the best probability so far, update the best edge info

  double probability = QInsideComputer::compute(a, dtr, log_normalisation_factor);

  // if(std::isnan(probability))
  //   return;

  if (probability > derivations[0].probability)
    {
      derivations[0].probability = probability;
      derivations[0].dtrs = &dtr;
    }

  write_scores(dtr, probability);
}


void MaxRuleProbabilityMultiple::finalize()
{
}


void MaxRuleProbabilityMultiple::pick_best_lexical(const LexicalDaughter & dtr)
{
  if(occ[&dtr] == nb_grammars) {
    packed_edge_probability_with_index p;
    p.dtrs = &dtr;

    const std::vector<AnnotationInfo>& upannots = get_annotations_backup();


    //    std::cout << "size annots: " << upannots.size() << std::endl;


    const LexicalDaughter* d = static_cast<const LexicalDaughter*>(p.dtrs);
    //          std::cout << *( d->get_rule()->get_coarser(upannots.size() - i + 1)) << std::endl;

    p.probability = 0;
    for (unsigned i = 0; i < upannots.size(); ++i)
      {
        if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
          continue;
        const std::vector<double>& rule_probs = d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

        double probability = 0;
        for(unsigned j = 0; j < rule_probs.size(); ++j) {
          if(!upannots[i].valid_prob_at(j, LorgConstants::NullProba)) continue;
          // std::cout << "out: " << upannots[i].outside_probabilities.array[j] << std::endl;
          // std::cout << "prb: " << rule_probs[j] << std::endl;
          probability += rule_probs[j] * upannots[i].outside_probabilities.array[j];
        }

        //        std::cout << probability << " " << std::log(probability) << " " << get_log_normalisation_factor(i) << std::endl;
        double res =
          // (get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity()) ?
          // -std::numeric_limits<double>::infinity() :
          std::log(probability) - get_log_normalisation_factor(i);
        // FIXME !!! We shouldn't resort to this hack
        //res = res > 0 ? 0 : res;
        p.probability += res;

        // if(p.probability > 0) {
        //   std::cout << *(d->get_rule()->get_coarser(upannots.size() - i - 1)) << std::endl;
        //   std::cout << p.probability << std::endl;

        //   std::cout << "TEST" << std::endl;
        //   std::cout << *(d->get_rule()) << std::endl;
        //   std::cout << *(d->get_rule()->get_coarser(0)) << std::endl;
        //   std::cout << *(d->get_rule()->get_coarser(1)) << std::endl;
        //   std::cout << *(d->get_rule()->get_coarser(2)) << std::endl;

        // }

      }

    // if(!(p.probability <= 0)) {
    //   std::cout << "pbl: " << p.probability << std::endl;
    // }

    assert(p.probability <=0.0001);
    if(p.probability > 0) p.probability = 0;


    if (candidates.empty() || p.probability > derivations[0].probability)
      derivations[0] = p;
    candidates.push_back(p);
  }

}


void MaxRuleProbabilityMultiple::pick_best_binary(const BinaryDaughter& dtr)
{

  if(occ[&dtr] == nb_grammars) {

    packed_edge_probability_with_index p;
    p.dtrs = &dtr;

    const std::vector<AnnotationInfo>& upannots = get_annotations_backup();





    //    std::cout << "binary case" << std::endl;

    const BinaryDaughter * d = static_cast<const BinaryDaughter*>(p.dtrs);

    PEdge& left  = d->left_pdaughter();
    const std::vector<AnnotationInfo>& leftannots = left.get_best().get_annotations_backup();

    PEdge& right = d->right_pdaughter();
    const std::vector<AnnotationInfo>& rightannots = right.get_best().get_annotations_backup();

    p.probability = 0;

    for (unsigned i = 0; i < upannots.size(); ++i)
      {
        if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
          continue;
        //          std::cout << *( d->get_rule()->get_coarser(upannots.size() - i + 1)) << std::endl;

        const std::vector<std::vector<std::vector<double> > >& rule_probs =
          d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

        p.probability += QInsideComputer::compute_simple(upannots[i],
                                                                                   get_log_normalisation_factor(i),
                                                                                   leftannots[i],
                                                                                   rightannots[i],
                                                                                   rule_probs);
      }

    p.probability += left.get_best().get(0).probability + right.get_best().get(0).probability;

    // if(!(p.probability <= 0)) {
    //   std::cout << "pbb: " << p.probability << std::endl;
    // }
    assert(p.probability <= 0.00001);
    if (p.probability>0) p.probability = 0;
    //      std::cout << p.probability << std::endl;

    //assert(p.probability != -std::numeric_limits<double>::infinity());


    if (candidates.empty() || p.probability > derivations[0].probability)
      derivations[0] = p;
    candidates.push_back(p);

  }
}


void MaxRuleProbabilityMultiple::pick_best_unary(const UnaryDaughter & dtr)
{
  if(occ[&dtr] == nb_grammars) {

    packed_edge_probability_with_index p;
    p.dtrs = &dtr;

    const std::vector<AnnotationInfo>& upannots = get_annotations_backup();


    const UnaryDaughter* d = static_cast<const UnaryDaughter*>(p.dtrs);

    LBEdge& left  = d->lbdaughter();

    if(left.get_best().get(0).dtrs) {

      const std::vector<AnnotationInfo>& leftannots = left.get_best().get_annotations_backup();

      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
            continue;
          const std::vector<std::vector<double> >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                                                     get_log_normalisation_factor(i),
                                                                                     leftannots[i],
                                                                                     rule_probs);
        }


      //          std::cout << p.probability << " " << left.get_prob_model().get(0).probability << std::endl;
      //std::cout << *(d->get_rule()) << std::endl;
      p.probability +=  left.get_best().get(0).probability;
      //          std::cout << p.probability << std::endl;
      assert(p.probability <= 0.0000001); // rounding errors !

    // if(!(p.probability <= 0)) {
    //   std::cout << "pbu: " << p.probability << std::endl;
    // }

      if(p.probability>0) p.probability = 0;

      //assert(p.probability != -std::numeric_limits<double>::infinity());

    }
    if (candidates.empty() || p.probability > derivations[0].probability)
      derivations[0] = p;

    candidates.push_back(p);
  }
}

//read scores and pick best
void MaxRuleProbabilityMultiple::pick_best()
{
  // free memory
  score_map_type().swap(scores);
  occ_map_type().swap(occ);

  if(!candidates.empty()) {
    if(candidates.size() > size) {
      //      std::cout << candidates.size() << std::endl;
      std::nth_element(candidates.begin(),candidates.begin()+size,candidates.end(), std::greater<packed_edge_probability_with_index>());
      candidates.resize(size);

    }
    std::make_heap(candidates.begin(),candidates.end());

    std::pop_heap(candidates.begin(),candidates.end());
    derivations[0] = candidates.back();
    candidates.pop_back();
  }
  else
    derivations.resize(0);

}


void MaxRuleProbabilityMultiple::backup_annotations(const AnnotationInfo& annotations)
{
  // std::cout << "outside " ;
  //   for (int i = 0; i < annotations.outside_probabilities.array.size(); ++i)
  //     {
  //       std::cout << annotations.outside_probabilities.array[i] << " " ;
  //     }
  //   std::cout << std::endl;

  //   std::cout << "inside " ;
  //   for (int i = 0; i < annotations.inside_probabilities.array.size(); ++i)
  //     {
  //       std::cout << annotations.inside_probabilities.array[i] << " " ;
  //     }
  //   std::cout << std::endl;
  annotations_backup.push_back(annotations);
}


inline
std::vector<AnnotationInfo>&
MaxRuleProbabilityMultiple::get_annotations_backup() {return annotations_backup;}

inline
const std::vector<AnnotationInfo>& MaxRuleProbabilityMultiple::get_annotations_backup() const {return annotations_backup;}

///////////////////////////


void MaxRuleProbabilityMultiple::extend_derivation(PEdge* edge, unsigned i, bool licence_unaries)
{
  if(derivations.size() == i) {
    return;
  }

  if(derivations.size() > 0) {

    packed_edge_probability_with_index& last = derivations[derivations.size() -1];

    //    std::cout << "last.probability " << last.probability << std::endl;

    assert(last.probability <= 0);
    // if (!(last.probability <=0)) {
    //   std::cout << derivations.size() << std::endl;
    //   std::cout << last.probability << std::endl << std::endl;
    //   last.probability = 0;
    // }

    find_succ(edge,last,licence_unaries);
    //    std::cout << "after find_succ" << std::endl;
  }

  if (!candidates.empty()) {

    //get next element from the candidates and append it to derivations
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


void MaxRuleProbabilityMultiple::find_succ(PEdge* edge, packed_edge_probability_with_index& pep, bool licence_unaries)
{
  if(pep.dtrs->is_lexical())  {
    //std::cout << "find_suc lex" << std::endl;
    return;}
  // binary -> extend left and right daughters


  const std::vector<AnnotationInfo>& upannots = edge->get_best().get_annotations_backup();

  if(pep.dtrs->is_binary()) {

    //    std::cout << "binary case" << std::endl;

    const BinaryDaughter* d = static_cast<const BinaryDaughter*>(pep.dtrs);

    PEdge& left  = d->left_pdaughter();
    const std::vector<AnnotationInfo>& leftannots = left.get_best().get_annotations_backup();

    PEdge& right = d->right_pdaughter();
    const std::vector<AnnotationInfo>& rightannots = right.get_best().get_annotations_backup();

    //extend to the left
    //    std::cout << "bin extending on the left" << std::endl;
    unsigned nextleft = pep.left_index + 1;
    left.extend_derivation(nextleft+1,true);

    // we haven't reached the expected number of solutions
    if(nextleft < left.get_best().n_deriv()) {

      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
            continue;
          //          std::cout << *( d->get_rule()->get_coarser(upannots.size() - i + 1)) << std::endl;

          const std::vector<std::vector<std::vector<double> > >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                                                     get_log_normalisation_factor(i),
                                                                                     leftannots[i],
                                                                                     rightannots[i],
                                                                                     rule_probs);
        }

      p.probability += left.get_best().get(p.left_index).probability + right.get_best().get(p.right_index).probability;


      assert(p.probability <= 0);

      //      std::cout << p.probability << std::endl;

      // TODO : Find a proper way to remove duplicates !
      if (std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()) {
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }

    }

    //extend to the right
    unsigned nextright = pep.right_index + 1;

    right.extend_derivation(nextright+1,true);

    if(nextright < right.get_best().n_deriv()) {
      //      std::cout << "bin extending on the right" << std::endl;


      packed_edge_probability_with_index p(pep);
      p.right_index = nextright;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
            continue;
          const std::vector<std::vector<std::vector<double> > >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                                                     get_log_normalisation_factor(i),
                                                                                     leftannots[i],
                                                                                     rightannots[i],
                                                                                     rule_probs);
        }

      p.probability += left.get_best().get(p.left_index).probability + right.get_best().get(p.right_index).probability;



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

    //    std::cout << "unary case" << std::endl;

    const UnaryDaughter* d = static_cast<const UnaryDaughter*>(pep.dtrs);

    LBEdge& left  = d->lbdaughter();
    const std::vector<AnnotationInfo>& leftannots = left.get_best().get_annotations_backup();

    //        std::cout << * d->get_rule() << std::endl;


    //extend to the left
    unsigned nextleft = pep.left_index + 1;

    left.extend_derivation(nextleft+1, false);

    if(nextleft < left.get_best().n_deriv() ) {
      //        std::cout << "un extending" << std::endl;
      packed_edge_probability_with_index p(pep);
      p.left_index = nextleft;
      p.probability = 0;


      for (unsigned i = 0; i < upannots.size(); ++i)
        {
          if(get_log_normalisation_factor(i) == -std::numeric_limits<double>::infinity())
            continue;
          const std::vector<std::vector<double> >& rule_probs =
            d->get_rule()->get_coarser(upannots.size() - i - 1)->get_probability();

          p.probability += QInsideComputer::compute_simple(upannots[i],
                                                                                     get_log_normalisation_factor(i),
                                                                                     leftannots[i],
                                                                                     rule_probs);
        }

      p.probability +=  left.get_best().get(p.left_index).probability;



      assert(p.probability <= 0);

      //            std::cout << p.probability << std::endl;


      if(std::find_if(candidates.begin(), candidates.end(), test_helper(p)) == candidates.end()){
        candidates.push_back(p);
        push_heap(candidates.begin(), candidates.end());
      }
    }
  }
}




#endif /* _MAXRULEMULTIPLEPROBABILITY_CPP_ */