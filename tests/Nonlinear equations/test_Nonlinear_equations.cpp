#include <fstream>
#include <gtest/gtest.h>
#include "Vychmaty/Nonlinear equations/Nonlinear_equations.hpp"

TEST(Testing, Newton){
    std::cout << keplerSolver(0.8, M_PI / 4, 1e2, 1e-15);
}

double func_2(double const x){
    return std::tan(x);
}

double func(double const x){
    return x - std::atan(std::sqrt(1 - x * x));
}


TEST(Testing, SIMM){
    double x, y, x0 = 0, y0 = 0;

    x = SIM_solve<double>(func, 0.4, x0, 1e3, 1e-6);
    y = func_2(x);
    std::cout << x << " " << y << std::endl;    // выводит "0.649889 0.760029", что и является правдой (посторил в десмосе, точность в 6 знаков правильная)
                                                // аналогично и для второго корня
}
