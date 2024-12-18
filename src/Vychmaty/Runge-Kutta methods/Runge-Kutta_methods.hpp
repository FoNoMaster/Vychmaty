#ifndef VYCHMATY_RUNGE_KUTTA_METHODS_HPP
#define VYCHMATY_RUNGE_KUTTA_METHODS_HPP

#include <cmath>
#include <array>
#include <vector>
#include "Eigen/Dense"
#include <iostream>



template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg& initialState,
        const typename RHS::Argument& endTime,
        double step,
        const RHS& rhs){

    std::vector<typename RHS::StateAndArg> result {initialState};
    typename RHS::StateAndArg state = initialState;
    Eigen::Vector<double, RHS::dim> u, tmp;
    Eigen::Vector<Eigen::Vector<double, RHS::dim>, Table::cColumn.size()> k;


    for (typename RHS::Argument t = initialState.arg + step; t <= endTime; t += static_cast<typename RHS::Argument>(step)) {

        k[0] = rhs.calc(state);

        for(int i = 1; i < k.size(); i++){
            tmp.setZero();
            for(int j = 0; j < i; j++){
                tmp += k[j] * Table::table[i][j];
            }

            k[i] = rhs.calc({state.state + step * tmp, state.arg + step * Table::cColumn[i]});
        }

        u = k[0] * Table::bString[0];
        for(int i = 1; i < k.size(); i++){
            u += k[i] * Table::bString[i];
        }
        u *= step;
        u += state.state;

        state = {u, t};
        result.push_back(state);

        if (t + step >= endTime) {
            step = endTime - t;
            if(step == 0)
                break;
        }
    }
    return result;
}


template<typename T>
T dot(const std::vector<T>& a, const std::vector<T>& b){
    T res = 0;
    for(std::size_t i = 0; i < a.size(); i++)
        res += a[i] * b[i];
    return res;
}

template<typename T>
T norm(const std::vector<T>& v){
    return sqrt(dot(v, v));
}

struct StepControl{
    double minStep;
    double maxStep;
    double tolerance;
    double initialStep;
};


template<typename Table, typename RHS>  // таблица бутчера и класс правой части f
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg& initialState,
        const typename RHS::Argument& endTime,
        const StepControl& stepControl,
        const RHS& rhs){

    std::vector<typename RHS::StateAndArg> result {initialState};
    typename RHS::StateAndArg state = initialState;
    Eigen::Vector<double, RHS::dim> u1, u2, tmp;
    Eigen::Vector<Eigen::Vector<double, RHS::dim>, Table::cColumn.size()> k;

    double step = stepControl.initialStep, error;


    for (typename RHS::Argument t = initialState.arg + step; t <= endTime; t += static_cast<typename RHS::Argument>(step)) {


        k[0] = rhs.calc(state);

        for(std::size_t i = 1; i < k.size(); i++){
            tmp.setZero();
            for(std::size_t j = 0; j < i; j++){
                tmp += k[j] * Table::table[i][j];
            }

            k[i] = rhs.calc({state.state + step * tmp, state.arg + step * Table::cColumn[i]});
        }

        u1 = k[0] * Table::bString_1[0];
        u2 = k[0] * Table::bString_2[0];

        for(std::size_t i = 1; i < k.size(); i++){
            u1 += k[i] * Table::bString_1[i];
            u2 += k[i] * Table::bString_2[i];
        }

        u1 *= step;
        u1 += state.state;

        u2 *= step;
        u2 += state.state;

        error = (u1 - u2).norm();

        if(error > stepControl.tolerance){
            t -= step;
        }


        step = std::min(stepControl.maxStep, std::max(stepControl.minStep, 0.95 * step * std::pow(stepControl.tolerance / error, 1. / (Table::order))));


        if(error <= stepControl.tolerance) {

            if (t + step >= endTime) {
                step = endTime - t;
                if(step == 0){
                    break;
                }
            }
            state = {u2, t};
            result.push_back(state);
        }
    }
    return result;
}



#endif //VYCHMATY_RUNGE_KUTTA_METHODS_HPP
