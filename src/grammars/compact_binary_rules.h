// -*- mode: c++ -*-
#ifndef _COMPACT_BINARY_RULES_H_
#define _COMPACT_BINARY_RULES_H_

#include <vector>

namespace compact_binary_rules {
  
  // this structure stores rules sharing attribute rhs1
  // begin and end are precomputed for faster loops  
  template<typename info>
  struct vector_rhs1 {
    int rhs1;
    typedef std::vector<info, std::allocator<info> > data;
    typedef typename data::const_iterator const_iterator;
    data rules;
    const_iterator _begin;
    const_iterator _end;

    vector_rhs1() : rhs1(), rules(), _begin(), _end() {}
  };
  
  // this structure stores vector_rhs1 sharing attribute rhs0
  // begin and end are precomputed for faster loops
  template<typename info>
  struct vector_rhs0 {
    int rhs0;
    typedef std::vector< vector_rhs1<info>, std::allocator< vector_rhs1<info> > > data;
    typedef typename data::const_iterator const_iterator;
    data vrhs1;
    const_iterator _begin;
    const_iterator _end;

    vector_rhs0() : rhs0(), vrhs1(), _begin(), _end() {}
  };
  
  // this structure stores vector_rhs0
  // begin and end are precomputed for faster loops
  template<typename info>
  struct vector_brules {
    typedef std::vector< vector_rhs0<info>, std::allocator< vector_rhs0<info> > > data;
    typedef typename data::iterator iterator;
    typedef typename data::const_iterator const_iterator;
    data vrhs0;
    const_iterator _begin;
    const_iterator _end;
    
    vector_brules() : vrhs0(), _begin(), _end() {}

    template<class BinaryRule> 
    static vector_brules * convert(const std::vector<BinaryRule, std::allocator<BinaryRule> >&);

    inline const_iterator begin() const { return _begin; }
    inline const_iterator end() const { return _end; }
  };
}

namespace std {
  template<typename info>
  inline typename compact_binary_rules::vector_brules<const info *>::const_iterator begin(typename compact_binary_rules::vector_brules<const info *>* const & v) { return v->_begin; }
  template<typename info>
  inline typename compact_binary_rules::vector_brules<const info *>::const_iterator end(typename compact_binary_rules::vector_brules<const info *>* const & v) { return v->_end; }
  template<typename info>
  inline typename compact_binary_rules::vector_rhs0<const info *>::const_iterator begin(const typename compact_binary_rules::vector_rhs0<const info *> & v) { return v._begin; }
  template<typename info>
  inline typename compact_binary_rules::vector_rhs0<const info *>::const_iterator end(const typename compact_binary_rules::vector_rhs0<const info *> & v) { return v._end; }
  template<typename info>
  inline typename compact_binary_rules::vector_rhs1<const info *>::const_iterator begin(const typename compact_binary_rules::vector_rhs1<const info *> & v) { return v._begin; }
  template<typename info>
  inline typename compact_binary_rules::vector_rhs1<const info *>::const_iterator end(const typename compact_binary_rules::vector_rhs1<const info *> & v) { return v._end; }
};
#endif /* _COMPACT_BINARY_RULES_H_ */
