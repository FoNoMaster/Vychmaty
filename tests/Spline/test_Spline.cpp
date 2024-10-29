#include <gtest/gtest.h>
#include "Spline/Spline.hpp"
#include <cmath>
#include <fstream>


TEST(Spline, test1){

    std::ofstream fout("estestvennyi.txt");
    fout.clear();

    int N = 5;
    double end = 10, X;

    double max_error, error, max_x;

    while(N <= 160) {
        max_error = 0;
        std::vector<double> points(N), values(N);
        for (int j = 0; j < N; j++) {
            points[j] = j * end / (N - 1);
            values[j] = std::exp(points[j]);
        }
        CubicSpline cubicSpline(points, values, 0, 0);

        for(int i = 0; i < 1000; i++){
            error = std::abs(std::exp(i * end / 999) - cubicSpline.interpolate(i * end / 999));
            if(error > max_error){
                max_error = error;
                max_x = i * end / 999;
            }
        }

        fout << std::log(N) << " " << std::log(max_error) << std::endl;
        std::cout << N << " " << max_x << std::endl;
        N += 1;
    }

//        std::vector<double> points(N), values(N);
//        for (int j = 0; j < N; j++) {
//            points[j] = j * end / (N - 1);
//            values[j] = std::exp(points[j]);
//        }
//        CubicSpline cubicSpline(points, values, 0, 0);
//
//        for(int i = 0; i < 1000; i++){
//            error = std::abs(std::exp(i * end / 999) - cubicSpline.interpolate(i * end / 999));
//            if(error > max_error){
//                max_error = error;
//                max_x = i * end / 999;
//            }
//            fout << i * end / 999 << " " << std::log(error) << std::endl;
//        }
//
//        std::cout << N << " " << max_x << std::endl;


    fout.close();

}