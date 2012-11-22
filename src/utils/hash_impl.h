// -*- mode: c++ -*-
#ifndef _HASH_IMPL_H_
#define _HASH_IMPL_H_


namespace std
{
template <class U, class V>
class hash<pair<U,V> > {
 public:
  inline size_t operator()(const pair<U,V> & p) const {
    size_t seed = 0;
    seed ^= std::hash<U>()(p.first) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    seed ^= std::hash<V>()(p.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);

    return seed;
  }
};
};


#endif /* _HASH_IMPL_H_ */
