#include <gtest/gtest.h>
#include "Numerical_Differetiation/Numerical_differentiation.hpp"
#include <cmath>
#include <fstream>


TEST(Diff, test1){
    std::array<double, 2> points {-1, 1};
    DerivativeCoef<double, 2> result = calcDerivativeCoef<double, 2, 1>(points);

    std::array<double, 2> solution {-0.5, 0.5};

    ASSERT_NEAR(result.centralCoef_, 0, 1e-15);

    for(int i = 0; i < points.size(); i++){
        ASSERT_NEAR(result.otherCoefs_[i], solution[i], 1e-15);
    }
}

TEST(Diff, test2){
    std::array<double, 2> points {1, 2};
    DerivativeCoef<double, 2> result = calcDerivativeCoef<double, 2, 1>(points);

    std::array<double, 2> solution {2, -0.5};

    ASSERT_NEAR(result.centralCoef_, -1.5, 1e-15);

    for(int i = 0; i < points.size(); i++){
        ASSERT_NEAR(result.otherCoefs_[i], solution[i], 1e-15);
    }
}

TEST(Diff, test3){
    std::array<double, 2> points {-1, 1};
    DerivativeCoef<double, 2> result = calcDerivativeCoef<double, 2, 2>(points);

    std::array<double, 2> solution {1, 1};

    ASSERT_NEAR(result.centralCoef_, -2, 1e-15);

    for(int i = 0; i < points.size(); i++){
        ASSERT_NEAR(result.otherCoefs_[i], solution[i], 1e-15);
    }
}

TEST(Diff_graphs, N3_5){
    std::ofstream fout("N5.txt");
    fout.clear();

    double x0 = 1, h, error, func_res;
    const unsigned int N = 5;
    std::array<double, N> points {-2, -1, 1, 2, 937469387};
    DerivativeCoef<double, N> result = calcDerivativeCoef<double, N, 2>(points);

    for(int i = 0; i < 16; i++){
        h = std::pow(10, -i);
        func_res = result.centralCoef_ * std::exp(x0);
        for(int j = 0; j < result.otherCoefs_.size(); j++){
            func_res += result.otherCoefs_[j] * std::exp(x0 + points[j] * h);
        }
        func_res = func_res / h / h;

        std::cout << h << std::endl;
        std::cout << std::exp(x0) << std::endl;
        std::cout << func_res << std::endl << std::endl;

        error = std::abs(std::exp(x0) - func_res);

        fout << -i << " " << std::log10(error) << std::endl;
    }

    fout.close();
}
