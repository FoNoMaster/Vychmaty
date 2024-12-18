#ifndef VYCHMATY_BACKWARD_DIFF_HPP
#define VYCHMATY_BACKWARD_DIFF_HPP


#include <array>
#include <vector>
#include "Vychmaty/Runge-Kutta methods/Runge-Kutta_methods.hpp"

struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size> alpha = {11. / 6., -3, 1.5, -1. / 3};
    static constexpr std::array<double, size> betta = {1. / 3., 1.5, -1, -1. / 6};
};


struct IntegrationParameters{
    double step;  // шаг интегрирования
    double epsilon;  // точность решения нелинейного уравнения
    std::size_t maxIter;  // максимальное количество итераций для решения нелинейного уравнения
};


template<typename BDF, typename RKTable, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg& initialState,
        const typename RHS::Argument& endTime,
        const IntegrationParameters& parameters,
        const RHS& rhs){
    std::vector<typename RHS::StateAndArg> result = integrate<RKTable>(initialState, initialState.arg + (BDF::size - 2) * parameters.step,
                                                                       parameters.step, rhs);

    typename RHS::StateAndArg state = result.back();
    Eigen::Vector<double, RHS::dim> u;


    for (typename RHS::Argument t = state.arg + parameters.step; t <= endTime - parameters.step; t += static_cast<typename RHS::Argument>(parameters.step)) {

        u = parameters.step * rhs.calc(state);
        for(std::size_t i = 1; i < BDF::size; i++){
            u -= result[result.size() - i].state * BDF::betta[i];
        }
        u /= BDF::betta[0];

        for(std::size_t i = 0; i < parameters.maxIter; i++){
            const Eigen::Vector<double, RHS::dim> u_tmp = u;

            u = parameters.step * rhs.calc({u, t});

            for(std::size_t j = 1; j < BDF::size; j++){
                u -= result[result.size() - j].state * BDF::alpha[j];
            }
            u /= BDF::alpha[0];

            if((u - u_tmp).norm() < parameters.epsilon){
                break;
            }
        }

        state = {u, t};
        result.push_back(state);

//        if (t + step >= endTime) {
//            step = endTime - t;
//            if(step == 0)
//                break;
//        }
    }
    return result;
}





#endif //VYCHMATY_BACKWARD_DIFF_HPP
