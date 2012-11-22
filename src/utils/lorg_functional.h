// -*- mode: c++ -*-
#ifndef _LORG_FUNCTIONAL_H_
#define _LORG_FUNCTIONAL_H_

#include <functional>

// transform a method into a function
template<class Class, class Return, class... Args>
std::function<Return(Class&, Args&... args)> toFunc( Return(Class::*f)(Args&... args) )  {
  return std::function<Return(Class&, Args&... args)> (f);
}

// transform a const method into a function
template<class Class, class Return, class... Args>
auto toFunc( Return(Class::*f)(Args&... args) const) -> std::function<Return(const Class&, Args&... args)> {
  return std::function<Return(const Class&, Args&... args)> (f);
}

// // transform a const method into a function
// template<class Class, class Return, class... Args>
// auto toFunc( Return(Class::*f)(const Args... args) const) -> std::function<Return(const Class&, const Args... args)> {
//   return std::function<Return(const Class&, const Args... args)> (f);
// }

// do nothing
template<class Return, class... Args>
auto toFunc( std::function<Return(Args&... args)> & f) -> std::function<Return(Args&... args)> & {
  return f;
}

#endif /* _LORG_FUNCTIONAL_H_ */
