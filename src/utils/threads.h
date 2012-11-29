#ifndef _UTILS_THREADS_H_
#define _UTILS_THREADS_H_

#ifdef USE_THREADS

#include <tbb/atomic.h>
namespace {
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
};

#endif

#endif
