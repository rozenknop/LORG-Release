// -*- mode: c++ -*-
#ifndef _LORG_FUNCTIONAL_H_
#define _LORG_FUNCTIONAL_H_

#include <functional>

// transform a method into a function
template<class Class, class Return, class... Args>
std::function<Return(Class&, Args&... args)> toFunc( Return(Class::*f)(Args&... args) )  {
    return std::function<Return(Class&, Args&... args)> (f);
}

#endif /* _LORG_FUNCTIONAL_H_ */
