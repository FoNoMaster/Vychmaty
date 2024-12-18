#include "Vychmaty/Multi-step_Methods/Backward_Diff.hpp"
#include <gtest/gtest.h>
#include <fstream>

struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = {{{0, 0, 0, 0},
                                                                              {0.5, 0, 0, 0},
                                                                              {0, 0.5, 0, 0},
                                                                              {0, 0, 1, 0}}};
    static constexpr std::array<double, stages> cColumn = {0, 0.5, 0.5, 1};
    static constexpr std::array<double, stages> bString = {1./6, 1./3, 1./3, 1./6};
};


class Func_1 {

public:

    static constexpr unsigned int dim = 1;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    [[nodiscard]] static Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) {
        return Eigen::Vector<double, dim>{stateAndArg.arg * stateAndArg.arg * stateAndArg.arg};
    }
};

double Func_1_res(const double& t){
    return t * t * t * t / 4;
}


TEST(Testing, BDF4){
    std::ofstream fout("BDF4.txt");
    fout.clear();

    Func_1 func_1;
    Eigen::Vector<double, 1> shish;
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Func_1::StateAndArg> res = integrate<BDF4, RK4Table>({{shish}, 0}, 5, {step, 1e-2, 10000}, func_1);
        for(auto & re : res){
            if(std::abs(Func_1_res(re.arg) - re.state[0]) > error){
                error = std::abs(Func_1_res(re.arg) - re.state[0]);
            }
        }
        std::cout << res.size() << " ";
        fout << std::log10(step) << " " << std::log10(error) << std::endl;
        step *= 0.6;
    }

    fout.close();
}

class Oscillator {

public:

    static constexpr unsigned int dim = 2;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    [[nodiscard]] static Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) {
        return Eigen::Vector<double, dim>{-stateAndArg.state[1], stateAndArg.state[0]};
    }
};


TEST(Testing, BDF4_2){
    std::ofstream fout("BDF4_2.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Oscillator::StateAndArg> res = integrate<BDF4, RK4Table>({{shish}, 0}, 5, {step, 1e-11, 10000}, func_2);
        for(auto & re : res){
            double dx = std::cos(re.arg) - re.state[0], dy = std::sin(re.arg) - re.state[1];
            if(std::sqrt(dx * dx + dy * dy) > error){
                error = std::sqrt(dx * dx + dy * dy);
            }
        }

        fout << std::log10(step) << " " << std::log10(error) << std::endl;
        step *= 0.6;
    }

    fout.close();
}

struct CF_BDF5 {
    static constexpr std::size_t size = 5;
    static constexpr std::array<double, size> betta{1. / 4, 2.5 / 3., 1.5, -0.5, 1. / 12.};
    static constexpr std::array<double, size> alpha{25. / 12., 4., -3., 4. / 3., -1. / 4.};
};

TEST(Testing, BDF4_3){
    std::ofstream fout("BDF5.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Oscillator::StateAndArg> res = integrate<CF_BDF5, RK4Table>({{shish}, 0}, 5, {step, 1e-10, 10000}, func_2);
        for(auto & re : res){
            double dx = std::cos(re.arg) - re.state[0], dy = std::sin(re.arg) - re.state[1];
            if(std::sqrt(dx * dx + dy * dy) > error){
                error = std::sqrt(dx * dx + dy * dy);
            }
        }

        fout << std::log10(step) << " " << std::log10(error) << std::endl;
        step *= 0.6;
    }

    fout.close();
}


struct CF_BDF6 {
    static constexpr std::size_t size = 6;
    static constexpr std::array<double, size> alpha{-65. / 12., 10., -5., 5. / 3., -1. / 4., 5.};
    static constexpr std::array<double, size> betta{300. / 137., -300. / 137., 200. / 137., -75. / 137., 12. / 137., 60. / 137.};
};


TEST(Testing, BDF6){
    std::ofstream fout("BDF6.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Oscillator::StateAndArg> res = integrate<CF_BDF6, RK4Table>({{shish}, 0}, 5, {step, 1e-10, 10000}, func_2);
        for(auto & re : res){
            double dx = std::cos(re.arg) - re.state[0], dy = std::sin(re.arg) - re.state[1];
            if(std::sqrt(dx * dx + dy * dy) > error){
                error = std::sqrt(dx * dx + dy * dy);
            }
        }

        fout << std::log10(step) << " " << std::log10(error) << std::endl;
        step *= 0.6;
    }

    fout.close();
}


double vx(double x, double y, double z){
    return 398600441588889 * (-x + 1.5 * 0.001082636023 * (5 * x * z * z / (x*x+y*y+z*z) - x)) / pow(x*x+y*y+z*z, 1.5);
}

double vy(double x, double y, double z){
    return 398600441588889 * (-y + 1.5 * 0.001082636023 * (5 * y * z * z / (x*x+y*y+z*z) - y)) / pow(x*x+y*y+z*z, 1.5);
}

double vz(double x, double y, double z){
    return 398600441588889 * (-z + 1.5 * 0.001082636023 * (5 * z * z * z / (x*x+y*y+z*z) - 3 * z)) / pow(x*x+y*y+z*z, 1.5);
}

class Orbit {

public:

    static constexpr unsigned int dim = 6;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg {
        State state;
        Argument arg;
    };

    /* Вычисляет правую часть ДУ - функцию f*/
    Eigen::Vector<double, dim> calc(const StateAndArg &stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state[3], stateAndArg.state[4], stateAndArg.state[5],
                                          vx(stateAndArg.state[0], stateAndArg.state[1], stateAndArg.state[2]),
                                          vy(stateAndArg.state[0], stateAndArg.state[1], stateAndArg.state[2]),
                                          vz(stateAndArg.state[0], stateAndArg.state[1], stateAndArg.state[2])};
    }
};


TEST(Testing, Orbit){
    std::ofstream fout("Orbit_time.txt");
    fout.clear();

    Orbit orbit;
    Eigen::Vector<double, 6> shish{6800000, 0, 0, 0, std::sqrt(398600441588889 / (2 * 6800000)), std::sqrt(398600441588889 / (2 * 6800000))};
    double error, step = 100;

    std::vector<Orbit::StateAndArg> res = integrate<BDF4, RK4Table>({{shish}, 0}, 8600400, {step, 1e-15, 10000}, orbit);

    for(auto & re : res){
        fout << re.arg << " " << std::sqrt(re.state[0] * re.state[0] + re.state[1] * re.state[1] + re.state[2] * re.state[2]) << std::endl;
    }

    fout.close();
}

