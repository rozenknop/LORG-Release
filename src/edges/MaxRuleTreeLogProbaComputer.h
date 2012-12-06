// // -*- mode: c++ -*-
#ifndef _MAXRULEUPDATER_H_
#define _MAXRULEUPDATER_H_


#include "PCKYAllCell.h"

#include <numeric>
#include <math.h>

template<class ProbaModel>
class MaxRuleTreeLogProbaComputer
{
public:
  typedef typename ProbaModel::Edge Edge ;
  typedef typename Edge::LexicalDaughter LexicalDaughter;
  typedef typename Edge::UnaryDaughter UnaryDaughter;
  typedef typename Edge::BinaryDaughter BinaryDaughter;

  static double compute (const AnnotationInfo & up_annotations,
                                            const BinaryDaughter & dtr,
                                            double normalisation_factor,
                                            unsigned left_idx =0, unsigned right_idx = 0)
  {
    //    std::cout << left_idx << " : " << right_idx << std::endl;
    //    std::cout << *(dtr.get_rule()) << std::endl;
    
    const Edge & left  = dtr.left_daughter ();
    const Edge & right = dtr.right_daughter ();
    
    double probability = 0.0;
    
    const scaled_array & left_inside  = left.get_annotations ().inside_probabilities;
    const scaled_array & right_inside = right.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector<std::vector<std::vector<double >>> & rule_probs = dtr.get_rule ()->get_probability ();
    
    
    // std::cout << "up_outside " ;
    // for (int i = 0; i < up_outside.array.size(); ++i)
    //   {
      //     std::cout << up_outside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    // std::cout << "left_inside " ;
    // for (int i = 0; i < left_inside.array.size(); ++i)
    //   {
      //     std::cout << left_inside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    
    // std::cout << "right_inside " ;
    // for (int i = 0; i < right_inside.array.size(); ++i)
    //   {
      //     std::cout << right_inside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    unsigned size = rule_probs.size ();
    for (unsigned i = 0; i < size; ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      double temp = 0;
      const std::vector < std::vector < double >>&rule_probs_i =
      rule_probs[i];
      unsigned size_i = rule_probs_i.size ();
      for (unsigned j = 0; j < size_i; ++j)
      {
        if (!left.valid_prob_at (j))
          continue;
        double inner = 0;
        const std::vector < double >&rule_probs_ij = rule_probs_i[j];
        unsigned size_ij = rule_probs_ij.size ();
        for (unsigned k = 0; k < size_ij; ++k)
        {
          if (right.valid_prob_at (k))
            inner += rule_probs_ij[k] * right_inside.array[k];
        }
        temp += left_inside.array[j] * inner;
      }
      probability += up_outside.array[i] * temp;
    }
    
    
    double res = (std::log (probability) - normalisation_factor)
                 + left.get_prob_model().get (left_idx).probability
                 + right.get_prob_model().get (right_idx).probability;
    
    //assert(res <= 0);
    
    if (res > 0)
    {
      // for(unsigned i = 0; i < up_annotations.outside_probabilities.array.size(); ++i)
      //   //        if(up_annotations.outside_probabilities.array[i] != LorgConstants::NullProba)
      //   std ::cout << i << ":" << up_annotations.outside_probabilities.array[i] << " " ;
      // std::cout << std::endl;
      
      // std::cout << std::log(probability) << std::endl;
      // std::cout << normalisation_factor << std::endl;
      //        std::cout << res
      //<< " : " << *dtr.get_rule()
      //        << std::endl;
    }
    
    
    return (res > 0) ? 0 : res;
    
    
  }
  
  static void compute_best_indexes (const AnnotationInfo & up_annotations,
                                    const BinaryDaughter & dtrs,
                                    double normalisation_factor,
                                    unsigned &left_idx, unsigned &right_idx)
  {
    const Edge & left  = dtrs.left_daughter ()->get_edge (dtrs.get_rule ()->get_rhs0 ());
    const Edge & right = dtrs.right_daughter ()->get_edge (dtrs.get_rule ()->get_rhs1 ());
    
    const scaled_array & left_inside  = left.get_annotations ().inside_probabilities;
    const scaled_array & right_inside = right.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < std::vector < double >>>&rule_probs =
    dtrs.get_rule ()->get_probability ();
    
    double max_prob = -std::numeric_limits < double >::infinity ();
    double base = - normalisation_factor 
                  +  left.get_prob_model().get (0).probability
                  + right.get_prob_model().get (0).probability;
    
    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      const std::vector < std::vector < double >>&rule_probs_i =
      rule_probs[i];
      
      double contrib_i = std::log (up_outside.array[i]) + base;
      
      for (unsigned j = 0; j < rule_probs_i.size (); ++j)
      {
        if (!left.valid_prob_at (j))
          continue;
        const std::vector < double >&rule_probs_ij = rule_probs_i[j];
        
        double contrib_j = contrib_i + std::log (left_inside.array[j]);
        
        for (unsigned k = 0; k < rule_probs_ij.size (); ++k)
        {
          if (!right.valid_prob_at (k))
            continue;
          double contrib = std::log(right_inside.array[k]) + std::log(rule_probs_ij[k]) + contrib_j;
          
          if (contrib > max_prob)
          {
            left_idx = j;
            right_idx = k;
            max_prob = contrib;
          }
        }
      }
    }
  }
  
  static double compute (const AnnotationInfo & up_annotations,
                                            const UnaryDaughter & dtrs,
                                            double normalisation_factor,
                                            unsigned left_idx = 0)
  {
    double probability = 0;
    
    const Edge & left = dtrs.left_daughter();
    
    const scaled_array & left_inside = left.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < double >>&rule_probs = dtrs.get_rule ()->get_probability ();
    
    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      double inner (0.0);
      for (unsigned j = 0; j < rule_probs[i].size (); ++j)
      {
        if (!left.valid_prob_at (j))
          continue;
        inner += rule_probs[i][j] * left_inside.array[j];
      }
      probability += up_outside.array[i] * inner;
    }
    
    
    //FIXME: this should not happen because chart is clean ???
    // only relevant in kmax parsing
    if (      //probability != 0
      //     &&
      left.get_prob_model().n_deriv () != 0)
    {
        
      double res = (std::log(probability) - normalisation_factor) +
                   left.get_prob_model().get (left_idx).probability;
        // assert(res <= 0);
        
        if (res > 0)
        {
          // for(unsigned i = 0; i < up_annotations.outside_probabilities.array.size(); ++i)
          //   //        if(up_annotations.outside_probabilities.array[i] != LorgConstants::NullProba)
          //   std ::cout << i << ":" << up_annotations.outside_probabilities.array[i] << " " ;
          // std::cout << std::endl;
          
          // std::cout << std::log(probability) << std::endl;
          // std::cout << normalisation_factor << std::endl;
          //        std::cout << res
          //                  << " : " << *dtrs.get_rule()
          //                  << std::endl;
        }
        
        
        
        
        return (res > 0) ? 0 : res;
      }
      else
        return -std::numeric_limits < double >::infinity ();
  }
  
  static void compute_best_indexes (const AnnotationInfo & up_annotations,
                                    const UnaryDaughter & dtrs,
                                    double normalisation_factor,
                                    unsigned &left_idx)
  {
    const Edge & left =
    dtrs.left_daughter ()->get_edge (dtrs.get_rule ()->get_rhs0 ());
    
    const scaled_array & left_inside = left.get_annotations ().inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    const std::vector < std::vector < double >>&rule_probs = dtrs.get_rule ()->get_probability ();
    
    double max_prob = -std::numeric_limits < double >::infinity ();
    
    double base = - normalisation_factor 
                  + left.get_prob_model().get (0).probability;
    
    
    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      const std::vector < double >&rule_probs_i = rule_probs[i];
      double contrib_i = std::log (up_outside.array[i]) + base;
      for (unsigned j = 0; j < rule_probs_i.size (); ++j)
      {
        if (!left.valid_prob_at (j))
          continue;
        double contrib =
        std::log (left_inside.array[j]) + std::log (rule_probs_i[j]) +
        contrib_i;
        if (contrib > max_prob)
        {
          left_idx = j;
          max_prob = contrib;
        }
      }
    }
  }
  
  static double compute (const AnnotationInfo & up_annotations,
                                            const LexicalRuleC2f * rule_ptr,
                                            double normalisation_factor)
  {
    double probability = 0.0;
    
    for (unsigned i = 0; i < rule_ptr->get_probability ().size (); ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      //std::cout << up_annotations.outside_probabilities.array[i] << std::endl;
        probability +=
        rule_ptr->get_probability ()[i] *
        up_annotations.outside_probabilities.array[i];
    }
    
    double res = std::log (probability) - normalisation_factor;
    
    // if(res < 0) {
      //   std::cout << res << " : " << *rule_ptr  << std::endl;
    // }
    
    
    // if(res > 0) {
      //   // for(unsigned i = 0; i < up_annotations.outside_probabilities.array.size(); ++i)
    //   //   //        if(up_annotations.outside_probabilities.array[i] != LorgConstants::NullProba)
    //   //     std ::cout << i << ":" << up_annotations.outside_probabilities.array[i] << " " ;
    //   // std::cout << std::endl;
    
    //   // std::cout << std::log(probability) << std::endl;
    //   // std::cout << normalisation_factor << std::endl;
    //   std::cout << res
    //     // << " : " << *rule_ptr
    //             << std::endl;
    // }
    // //    assert(res <= 0);
    
    return res > 0 ? 0 : res;
  }
  
  
  static  double compute_simple (const AnnotationInfo & up_annotations,
                                            double normalisation_factor,
                                            const AnnotationInfo &
                                            left_annotations,
                                            const AnnotationInfo &
                                            right_annotations,
                                            const std::vector<std::vector<std::vector<double>>>&rule_probs)
  {
    double probability = 0.0;
    
    const scaled_array & left_inside = left_annotations.inside_probabilities;
    const scaled_array & right_inside =
    right_annotations.inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    
    // std::cout << "up_outside " ;
    // for (int i = 0; i < up_outside.array.size(); ++i)
    //   {
      //     std::cout << up_outside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    // std::cout << "left_inside " ;
    // for (int i = 0; i < left_inside.array.size(); ++i)
    //   {
      //     std::cout << left_inside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    
    // std::cout << "right_inside " ;
    // for (int i = 0; i < right_inside.array.size(); ++i)
    //   {
      //     std::cout << right_inside.array[i] << " " ;
    //   }
    // std::cout << std::endl;
    
    
    unsigned size = rule_probs.size ();
    for (unsigned i = 0; i < size; ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      double temp = 0;
      const std::vector < std::vector < double >>&rule_probs_i =
      rule_probs[i];
      unsigned size_i = rule_probs_i.size ();
      for (unsigned j = 0; j < size_i; ++j)
      {
        if (!left_annotations.valid_prob_at (j, LorgConstants::NullProba))
          continue;
        double inner = 0;
        const std::vector < double >&rule_probs_ij = rule_probs_i[j];
        unsigned size_ij = rule_probs_ij.size ();
        for (unsigned k = 0; k < size_ij; ++k)
        {
          if (right_annotations.
            valid_prob_at (k, LorgConstants::NullProba))
            inner += rule_probs_ij[k] * right_inside.array[k];
          //          std::cout << "inner " << inner << std::endl;
        }
        temp += left_inside.array[j] * inner;
      }
      probability += up_outside.array[i] * temp;
      //      std::cout << probability << std::endl;
    }
    
    double res = (std::log (probability) - normalisation_factor);
    
    //    std::cout << std::log(probability) << " " << normalisation_factor << std::endl;
    
    //    std::cout << "umps bin: " << res << std::endl;
    
    return (res > 0) ? 0 : res;
    
    
  }
  
  static double compute_simple (const AnnotationInfo & up_annotations,
                                            double normalisation_factor,
                                            const AnnotationInfo &
                                            left_annotations,
                                            const std::vector <
                                            std::vector <
                                            double > >&rule_probs)
  {
    double probability = 0;
    
    const scaled_array & left_inside = left_annotations.inside_probabilities;
    const scaled_array & up_outside = up_annotations.outside_probabilities;
    
    for (unsigned i = 0; i < rule_probs.size (); ++i)
    {
      if (!up_annotations.valid_prob_at (i, LorgConstants::NullProba))
        continue;
      double inner (0.0);
      for (unsigned j = 0; j < rule_probs[i].size (); ++j)
      {
        if (!left_annotations.valid_prob_at (j, LorgConstants::NullProba))
          continue;
        inner += rule_probs[i][j] * left_inside.array[j];
      }
      probability += up_outside.array[i] * inner;
    }
    
    
    //FIXME: this should not happen because chart is clean ???
    // only relevant in kmax parsing
    // if(//probability != 0
    //    //     &&
    //    left.get_prob_model().n_deriv() != 0
    //    ) {
      
      double res = (std::log (probability) - normalisation_factor);
      return (res > 0) ? 0 : res;
      // }
      // else
      //   return - std::numeric_limits<double>::infinity();
  }
};

#endif /* _MAXRULEPROBABILITY_H_ */
