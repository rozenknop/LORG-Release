#ifndef lorg_functional_h
#define lorg_functional_h

#include <functional>

template<class Class, class Return, class... Args>
auto myMethod( Class * object, Return(Class::*f)(Args&&... args)) -> function<Return(Args&&... args)> {
    return [object,f](Args&&... args){return (object->*f)(std::forward<Args&&>(args)...);};
}

template<class Class, class Return, class... Args>
std::function<Return(Class&, Args&... args)> toFunc( Return(Class::*f)(Args&... args) )  {
    return std::function<Return(Class&, Args&... args)> (f);
}

 
#endif