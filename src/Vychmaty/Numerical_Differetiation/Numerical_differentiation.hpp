#ifndef VYCHMATY_NUMERICAL_DIFFERENTIATION_HPP
#define VYCHMATY_NUMERICAL_DIFFERENTIATION_HPP

#include "Eigen/Dense"
#include <array>


template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef_;
    std::array<RealType, N> otherCoefs_;

    DerivativeCoef(const RealType& centralCoef, const std::array<RealType, N>& otherCoefs):centralCoef_(centralCoef), otherCoefs_(otherCoefs){}

};

template<typename RealType, unsigned int N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {

    Eigen::MatrixX<RealType> A(N + 1, N + 1);
    A.setZero();
    A(0, 0) = 1;

    Eigen::VectorX<RealType> b(N + 1);
    b.setZero();
    b(L) = 1;


    for(unsigned int j = 1; j < N + 1; j++){
        A(0, j) = 1;
        for(unsigned int i = 1; i < N + 1; i++){
            A(i, j) = A(i - 1, j) * points[j - 1] / i;
        }
    }

    Eigen::VectorX<RealType> x = A.lu().solve(b);

    std::cout << A << std::endl;
//    std::cout << b << std::endl << std::endl;
    std::cout << x << std::endl << std::endl << std::endl;

    std::array<RealType, N> otherCoefs;
    for(int i = 1; i < N + 1; i++){
        otherCoefs[i - 1] = x(i);
    }

    return DerivativeCoef<RealType, N> (x(0), otherCoefs);
}


#endif //VYCHMATY_NUMERICAL_DIFFERENTIATION_HPP
