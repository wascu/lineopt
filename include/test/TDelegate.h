//
// Created by wasku on 19-7-27.
//

#ifndef EIGENSAMPLE_TDELEGATE_H
#define EIGENSAMPLE_TDELEGATE_H

template <class R, typename... Args>
class  TDelegate
{
public:
    TDelegate( R  (*f)(Args...) ):m_f(f) {}

    R operator()(Args&&... args)
    {
        return (*m_f)(std::forward<Args>(args) ...);
    }

private:
    R  (*m_f)(Args...);
};

template <class T, class R, typename... Args>
TDelegate<T, R, Args...> CreateDelegate(T* t, R (T::*f)(Args...))
{
    return TDelegate<T, R, Args...>(t, f);
}

template <class R, typename... Args>
TDelegate< R, Args...> CreateDelegate( R (*f)(Args...))
{
    return TDelegate<R, Args...>( f);
}

class TTW{
public:
        void test();
    };

#endif //EIGENSAMPLE_TDELEGATE_H
