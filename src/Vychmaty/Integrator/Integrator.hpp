#ifndef VYCHMATY_INTEGRATOR_HPP
#define VYCHMATY_INTEGRATOR_HPP


#include <iostream>
#include <array>
#include <type_traits>
#include <cmath>


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());


template <typename Callable>
struct Node {
    typename ArgumentGetter<Callable>::Argument point;
    double weight;
};

template <typename Callable, std::size_t N>
struct Nodes {
    static std::array<Node<Callable>, N> nodes;
};

template <typename Callable>
struct Nodes<Callable, 3> {
    static constexpr std::array<Node<Callable>, 3> nodes = {
            Node<Callable>{-0.77459666924148338, 0.5555555555555556},
            Node<Callable>{0.00000000000000000, 0.8888888888888889},
            Node<Callable>{0.77459666924148338, 0.5555555555555556}
    };
};

template <typename Callable>
struct Nodes<Callable, 4> {
    static constexpr std::array<Node<Callable>, 4> nodes = {
                    Node<Callable>{-0.8611363115940526, 0.3478548451374539},
                    Node<Callable>{-0.3399810435848563, 0.6521451548625461},
                    Node<Callable>{0.3399810435848563, 0.6521451548625461},
                    Node<Callable>{0.8611363115940526, 0.3478548451374539}
            };
};

template <typename Callable>
struct Nodes<Callable, 8> {
    static constexpr std::array<Node<Callable>, 8> nodes = {
                    Node<Callable>{-0.9602898564975362, 0.10122853629038},
                    Node<Callable>{-0.7966664774136267, 0.22238103445337},
                    Node<Callable>{-0.5255324099163290, 0.31370664587789},
                    Node<Callable>{-0.1834346424956498, 0.36268378337836},
                    Node<Callable>{0.1834346424956498, 0.36268378337836},
                    Node<Callable>{0.5255324099163290, 0.31370664587789},
                    Node<Callable>{0.7966664774136267, 0.22238103445337},
                    Node<Callable>{0.9602898564975362, 0.10122853629038}
            };
};

/* Функция производит интегрирование на одном отрезке */
template<std::size_t N, typename Callable>
decltype(auto) integrate(
        const Callable& func,  // Интегрируемая функция
        const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
        const typename ArgumentGetter<Callable>::Argument& end  // конец отрезка
){
    typename ArgumentGetter<Callable>::Argument sum = 0;
    for (int i = 0; i < N; i++) {
        sum += Nodes<Callable, N>::nodes[i].weight * func((start + end) / 2 + (end - start) / 2 * Nodes<Callable, N>::nodes[i].point);
    }
    return sum * (end - start) / 2;
}

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<std::size_t N, typename Callable>
decltype(auto) integrate(
        const Callable& func,  // Интегрируемая функция
        const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
        const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
        const Dif<typename ArgumentGetter<Callable>::Argument>& dx  // Длина подотрезка
){
    typename ArgumentGetter<Callable>::Argument I = 0;
    for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dx)) {
        if (x + dx >= end) {
            I += integrate<N>(func, x, end);
            break;
        }
        I += integrate<N>(func, x, x + dx);
    }
    return I;
}

template<std::size_t N, typename Callable>
decltype(auto) integrateFast(
        const Callable& func,  // Интегрируемая функция
        const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
        const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
        const typename ArgumentGetter<Callable>::Argument& precision  // Точность
){
    typename ArgumentGetter<Callable>::Argument error;
    typename ArgumentGetter<Callable>::Argument I1 = 0, I2 = 0;
    Dif<typename ArgumentGetter<Callable>::Argument> dx = end - start;

    for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dx)) {
        I1 += integrate<N>(func, x, x + dx);
    }

    dx /= 2;

    for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dx)) {
        I2 += integrate<N>(func, x, x + dx);
    }

    error = (I1 - I2) / (std::pow(2, 2 * N) - 1);

    std::cout << error << std::endl;

    while (std::abs(error) > precision){
        I1 = I2;
        I2 = 0;

        dx /= 2;

        for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dx)) {
            I2 += integrate<N>(func, x, x + dx);
        }

        error = (I1 - I2) / (std::pow(2, 2 * N) - 1);

        std::cout << error << std::endl;
    }
    return I2 - error;
}


template<std::size_t N, typename Callable>
decltype(auto) integrate(
        const Callable& func,  // Интегрируемая функция
        const typename ArgumentGetter<Callable>::Argument& start,  // начало отрезка
        const typename ArgumentGetter<Callable>::Argument& end,  // конец отрезка
        const Dif<typename ArgumentGetter<Callable>::Argument>& dx,  // Длина подотрезка
        const typename ArgumentGetter<Callable>::Argument& precision,  // Точность
        const unsigned int& iterations // Максимальное количество итераций
){
    typename ArgumentGetter<Callable>::Argument error;
    typename ArgumentGetter<Callable>::Argument I1 = 0, I2 = 0;
    Dif<typename ArgumentGetter<Callable>::Argument> dxx = dx;

    for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dxx)) {
        if (x + dxx >= end) {
            I1 += integrate<N>(func, x, end);
            break;
        }
        I1 += integrate<N>(func, x, x + dxx);
    }

    dxx /= 2;

    for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dxx)) {
        if (x + dxx >= end) {
            I2 += integrate<N>(func, x, end);
            break;
        }
        I2 += integrate<N>(func, x, x + dxx);
    }

    error = (I1 - I2) / (std::pow(2, 2 * N) - 1);

    unsigned int it = 0;
    while (std::abs(error) > precision && it < iterations){
        I1 = I2;
        I2 = 0;

        dxx /= 2;

        for (typename ArgumentGetter<Callable>::Argument x = start; x <= end; x += static_cast<typename ArgumentGetter<Callable>::Argument>(dxx)) {
            if (x + dxx >= end) {
                I2 += integrate<N>(func, x, end);
                break;
            }
            I2 += integrate<N>(func, x, x + dxx);
        }

        error = (I1 - I2) / (std::pow(2, 2 * N) - 1);

        std::cout << error << std::endl;

        it++;
    }
    return I2 - error;
}


#endif //VYCHMATY_INTEGRATOR_HPP
