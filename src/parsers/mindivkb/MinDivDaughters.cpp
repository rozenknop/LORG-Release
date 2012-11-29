#include "edges/PackedEdge.h"
#include "ChartCKY.h"
#include "parsers/mindivkb/ParserCKYAllMinDivKB.h"

#include "MinDivDaughters.h"

#ifdef USE_THREADS

#include <tbb/atomic.h>
using tbb::atomic;

inline void operator+=(atomic<double>& un, const double & deux)
/**
 * @brief atomic += operator used for concurrent computations of outside probabilities
 *
 * @param un left value to be modified
 * @param deux additive term
 * @return void
 **/
{
  double oldx, newx ;
  do {
    oldx = un ;
    newx = oldx+deux ;
  } while (un.compare_and_swap(newx,oldx) != oldx);
}

#endif





double
MinDivBRule::update_outside_annotations_return_marginal(const std::vector< double >& up_out, 
                                                        const std::vector< double >& left_in, 
                                                        const std::vector< double >& right_in, 
                                                        std::vector< double >& left_out, 
                                                        std::vector< double >& right_out) const
{
  double marginal = 0.;
  #ifdef USE_THREADS
  std::vector<atomic<double>> & lo = *reinterpret_cast<std::vector<atomic<double>> *>(&left_out);
  std::vector<atomic<double>> & ro = *reinterpret_cast<std::vector<atomic<double>> *>(&right_out);
  #else
  std::vector<double> & lo = left_out ;
  std::vector<double> & ro = right_out ;
  #endif
  for(unsigned short i = 0; i < probabilities.size(); ++i) {
    if(up_out[i] == LorgConstants::NullProba || up_out[i] == 0.0) continue;
    const std::vector<std::vector<double> >& dim_i = probabilities[i];
    for(unsigned short j = 0; j < dim_i.size(); ++j) {
      const std::vector<double>& dim_j = dim_i[j];
      double temp4left = 0.0;
      double factor4right = 0.0;
      if(left_in[j] != LorgConstants::NullProba) factor4right = up_out[i] * left_in[j];
      for(unsigned short k = 0; k < dim_j.size(); ++k) {
        const double& t = dim_j[k];
        // if(right_in[k] != LorgConstants::NullProba) temp4left += right_in[k] * t;
        // if(right_out[k] != LorgConstants::NullProba) right_out[k] += factor4right * t;
        
        // I and O are always Null at the same time
        if(right_in[k] != LorgConstants::NullProba) {
          temp4left += right_in[k] * t;
          ro[k] += factor4right * t;
        }
      }
      if(lo[j] != LorgConstants::NullProba) {
        double delta_left = up_out[i] * temp4left;
        lo[j] += delta_left;
        marginal += delta_left * left_in[j];
      }
    }
  }
  return marginal ;
}


double MinDivURule::update_outside_annotations_return_marginal(const std::vector< double >& up, 
                                                               const std::vector< double >& in_left, 
                                                               std::vector< double >& out_left) 
const
{
  double marginal = 0.0;
  for(unsigned short i = 0 ; i < probabilities.size();++i) {
    if(up[i] == LorgConstants::NullProba) continue;
    const std::vector<double>& dim_i = probabilities[i];
    for(unsigned short j = 0 ; j < dim_i.size();++j) {
      if(out_left[j] == LorgConstants::NullProba) continue;
      double delta = up[i] * dim_i[j] ;
      out_left[j] += delta ;
      marginal += delta * in_left[j] ;
    }
  }
  return marginal;
}


double MinDivLRule::update_outside_annotations_return_marginal(const std::vector< double >& up) const
{
  double marginal = 0;
  for(unsigned i = 0 ; i < probabilities.size();++i) {
    if(up[i] == LorgConstants::NullProba) continue;
    //if( up[i] == 0 ||      probabilities[i] == 0) continue;
    double delta = up[i] * probabilities[i];
    marginal += delta;
  }
  return marginal;
}


