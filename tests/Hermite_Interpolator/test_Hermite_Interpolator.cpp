#include <gtest/gtest.h>
#include "Hermite_Interpolator/Hermite_Interpolator.hpp"
#include <cmath>
#include <fstream>


TEST(Hermite_Interpolation, N3_5){
    std::ofstream fout("N5.txt");
    fout.clear();

    const unsigned int N = 5;
    std::array<double, N> points, values, deriv;
    double end = 2, max_error, error;

    for(int k = 0; k <= 5; k++){
        max_error = 0;
        for(int i = 0; i < N; i++){
            points[i] = i * end / (N - 1);
            values[i] = std::exp(points[i]);
            deriv[i] = std::exp(points[i]);
        }
        HermiteInterpolator<double, double, N> interpolator(points, values, deriv);

        for(int i = 0; i < 1000; i++){
            error = std::abs(std::exp(i * end / 999) - interpolator.interpolate(i * end / 999));
            if(error > max_error){
                max_error = error;
            }
        }
        std::cout << max_error << std::endl;
        fout << std::log10(end) << " " << std::log10(max_error) << std::endl;
        end /= 2;
    }

    fout.close();
}
