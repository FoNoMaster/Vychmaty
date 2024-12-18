#ifndef VYCHMATY_NONLINEAR_EQUATIONS_HPP
#define VYCHMATY_NONLINEAR_EQUATIONS_HPP

#include <cmath>
#include <fstream>


double keplerSolver(const double ecc, const double meanAnomaly, const unsigned int maxIter, const double tol){
    double E = meanAnomaly;
    for(unsigned int i = 0; i < maxIter; i++){
        const double dE = - (E - ecc * std::sin(E) - meanAnomaly) / (1 - ecc * std::cos(E));
        E += dE;
        if(std::abs(dE) < tol){
            return E;
        }
    }
    return E;
}

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};


template<typename RealType, typename Callable>
decltype(auto) SIM_solve(
        const Callable& func,                                             // функция F
        const RealType& tau,                                              // шаг тау
        const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
        const unsigned int& nIteration                                    // количество итераций
){
    typename ArgumentGetter<Callable>::Argument x_new = initialGuess, x_old;
    for(unsigned int i = 0; i < nIteration; i++){
        x_old = x_new;
        x_new = x_old - tau * func(x_old);
    }
    return x_new;
}

template<typename RealType, typename Callable>
decltype(auto) SIM_solve(
        const Callable& func,                                             // функция F
        const RealType& tau,                                              // шаг тау
        const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
        const unsigned int& nIteration,                                     // количество итераций
        const RealType& tol
){
    typename ArgumentGetter<Callable>::Argument x_new = initialGuess, x_old;
    for(unsigned int i = 0; i < nIteration; i++){
        x_old = x_new;
        x_new = x_old - tau * func(x_old);

        if(std::abs(x_old - x_new) <= tol){
            return x_new;
        }
    }
    return x_new;
}



#endif //VYCHMATY_NONLINEAR_EQUATIONS_HPP
