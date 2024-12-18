#include "Vychmaty/Runge-Kutta methods/Runge-Kutta_methods.hpp"
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


struct  BT_DP45 {
    static constexpr std::size_t stages = 7;
    static constexpr std::size_t order = 5;
    static constexpr std::array<std::array<double, stages>, stages> table{{{0, 0, 0, 0, 0, 0, 0},
                                                                           {1. / 5., 0, 0, 0, 0, 0, 0},
                                                                           {3. / 40., 9. / 40., 0, 0, 0, 0, 0},
                                                                           {44. / 45., -56. / 15., 32. / 9., 0, 0, 0, 0},
                                                                           {19372. / 6561., -25360. / 2187., 64448. / 6561., -212. / 729., 0, 0, 0},
                                                                           {9017. / 3168., -355. / 33., 46732. / 5247., 49. / 176., -5103. / 18656., 0, 0},
                                                                           {35. / 384., 0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0}}};
    static constexpr std::array<double, stages> cColumn{0, 1. / 5., 3. / 10., 4. / 5., 8. / 9., 1., 1.};
    static constexpr std::array<double, stages> bString_1{35. / 384., 0, 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0};
    static constexpr std::array<double, stages> bString_2{5179. / 57600., 0, 7571. / 16695., 393. / 640., -92097. / 339200., 187. / 2100., 1. / 40.};
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


TEST(Testing, Runge_Kutta_1){
    std::ofstream fout("RK4.txt");
    fout.clear();

    Func_1 func_1;
    Eigen::Vector<double, 1> shish;
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Func_1::StateAndArg> res = integrate<RK4Table>({{shish}, 0}, 5, step, func_1);
        for(auto & re : res){
            if(std::abs(Func_1_res(re.arg) - re.state[0]) > error){
                error = std::abs(Func_1_res(re.arg) - re.state[0]);
            }
        }
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


TEST(Testing, Runge_Kutta_2){
    std::ofstream fout("RK4_2.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    double error, step = 1;

    for(int i = 0; i < 25; i++) {
        error = 0;
        std::vector<Oscillator::StateAndArg> res = integrate<RK4Table>({{shish}, 0}, 5, step, func_2);
        for(auto & re : res){
            double dx = std::cos(re.arg) - re.state[0], dy = std::sin(re.arg) - re.state[1];
            if(std::sqrt(dx * dx + dy * dy) > error){
                error = std::sqrt(dx * dx + dy * dy);
            }
        }

//        for(auto & re : res){
//            double dx = std::abs(std::cos(re.arg) - re.state[0]);
//            if(std::sqrt(dx * dx) > error){
//                error = std::sqrt(dx * dx);
//            }
//        }

//        for(auto & re : res){
//            if(re.state[0] * re.state[0] + re.state[1] * re.state[1] - 1 > error){
//                error = re.state[0] * re.state[0] + re.state[1] * re.state[1] - 1;
//            }
//        }

        fout << std::log10(step) << " " << std::log10(error) << std::endl;
        step *= 0.6;
    }

    fout.close();
}


TEST(Testing, Runge_Kutta_23){
    std::ofstream fout("RK4_2_3.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    double step = 3e-2;
    std::vector<Oscillator::StateAndArg> res = integrate<RK4Table>({{shish}, 0}, 5, step, func_2);

    for(auto & re : res) {
        fout.precision(10);
        fout << "time: " << re.arg << "   " << re.state[0] << " " << re.state[1] << "  " << "real solution: " << std::cos(re.arg) << " " << std::sin(re.arg)
         << " delta: " << std::sqrt((std::cos(re.arg) - re.state[0]) * (std::cos(re.arg) - re.state[0]) + (std::sin(re.arg) - re.state[1]) * (std::sin(re.arg) - re.state[1])) << std::endl;
    }

    fout.close();
}

//TEST(Testing, DP_23){
//    std::ofstream fout("DP4_2_3.txt");
//    fout.clear();
//
//    Oscillator func_2;
//    Eigen::Vector<double, 2> shish{1, 0};
//    double step = 1e-4;
//    std::vector<Oscillator::StateAndArg> res = integrate<BT_DP45>({{shish}, 0}, 5, step, func_2);
//
//    for(auto & re : res) {
//        fout.precision(15);
//        fout << re.arg << " " << re.state[0] << " " << re.state[1] << std::endl;
//    }
//
//    fout.close();
//}


double R1(const double& x, const double& y, const double& m){
    return std::pow((x + m) * (x + m) + y * y, 1.5);
}

double R2(const double& x, const double& y, const double& M){
    return std::pow((x - M) * (x - M) + y * y, 1.5);
}

double u(const double& v, const double& x, const double& y, const double& m, const double& M){
    return x + 2 * v - M * (x + m) / R1(x, y, m) - m * (x - M) / R2(x, y, M);
}

double v(const double& u, const double& x, const double& y, const double& m, const double& M){
    return y - 2 * u - y * (M / R1(x, y, m) + m / R2(x, y, M));
}


class Arenstorf {

public:

    static constexpr unsigned int dim = 4;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    [[nodiscard]] static Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) {
        return Eigen::Vector<double, dim>{u(stateAndArg.state[1], stateAndArg.state[2], stateAndArg.state[3], 0.012277471, 0.987722529),
                                          v(stateAndArg.state[0], stateAndArg.state[2], stateAndArg.state[3], 0.012277471, 0.987722529),
                                          stateAndArg.state[0],
                                          stateAndArg.state[1]};
    }
};



TEST(Testing, Arenstorf){
    std::ofstream fout("BD.txt");
    fout.clear();
    fout.precision(20);

    Arenstorf func_3;
    Eigen::Vector<double, 4> shish{0, -2.031732629557337, 0.994, 0};
    StepControl stepControl{1e-6, 1, 1e-15, 1e-2};
    std::vector<Arenstorf::StateAndArg> res = integrate<BT_DP45>({{shish}, 0}, 11.124340337, stepControl, func_3);

    for(auto & re : res) {
        fout << re.state[2] << " " << re.state[3] << std::endl;
    }

    fout.close();
}

TEST(Testing, Arenstorf_RK){
    std::ofstream fout("RK_Arenstorf.txt");
    fout.clear();

    Arenstorf func_3;
    Eigen::Vector<double, 4> shish{0, -2.031732629557337, 0.994, 0};
    double step = 1e-3;
    std::vector<Arenstorf::StateAndArg> res = integrate<RK4Table>({{shish}, 0}, 5, step, func_3);

    for(auto & re : res) {
        fout << re.state[2] << " " << re.state[3] << std::endl;
    }

    fout.close();
}

TEST(Testing, BT_DP45){
    std::ofstream fout("BT_DP45.txt");
    fout.clear();

    Oscillator func_2;
    Eigen::Vector<double, 2> shish{1, 0};
    StepControl stepControl{1e-8, 1, 1e-15, 1};
    std::vector<Oscillator::StateAndArg> res = integrate<BT_DP45>({{shish}, 0}, 5, stepControl, func_2);

    for(auto & re : res) {
        fout.precision(20);
        fout << re.state[0] << " " << re.state[1] << std::endl;
    }

    fout.close();
}

