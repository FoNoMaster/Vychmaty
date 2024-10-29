#include "Integrator/Integrator.hpp"
#include <gtest/gtest.h>
#include <fstream>

double sinn(double const x){
    return std::sin(x);
}

TEST(qdwqd, lprwlb){
    std::ofstream fout("integrator3.txt");
    fout.clear();

    double start = 0, end = 10;
    double error, dx;

    for(int i = 1; i <= 100; i++){
        dx = 10. / i;
        error = std::abs(-std::cos(10) + std::cos(0) - integrate<3>(sinn, start, end, dx));
        fout << std::log10(dx) << " " << std::log10(error) << std::endl;
    }

    fout.close();
}

TEST(orwas, sodfgk){
    std::ofstream fout("integrator8_epsilon.txt");
    fout.clear();

    double start = 0, end = 10;
    double error, precision;

    for(int i = 0; i <= 150; i++){
        precision = std::pow(10, -i / 10);
        error = std::abs(-std::cos(10) + std::cos(0) - integrate<8>(sinn, start, end, 10, precision, 10000));
        fout << std::log10(precision) << " " << std::log10(error) << std::endl;
    }

    fout.close();
}