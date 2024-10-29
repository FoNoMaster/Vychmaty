#ifndef VYCHMATY_SPLINE_HPP
#define VYCHMATY_SPLINE_HPP


#include <vector>
#include <type_traits>
#include <iostream>
#include "Eigen/Dense"


template<typename T>
class ThreeDiagonalMatrix{
private:
    std::vector<T> a_;
    std::vector<T> b_;
    std::vector<T> c_;
public:
    ThreeDiagonalMatrix(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c): a_(a), b_(b), c_(c){}

    const std::vector<T>& get_a() const{return a_;}
    const std::vector<T>& get_b() const{return b_;}
    const std::vector<T>& get_c() const{return c_;}

    const T& a(const std::size_t& i) const{return a_[i];}
    const T& b(const std::size_t& i) const{return b_[i];}
    const T& c(const std::size_t& i) const{return c_[i];}
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());


/** Функция для решения методм  прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> ThreeDiagonalSolve(const ThreeDiagonalMatrix<mType>& A, const std::vector<cType>& d){
    std::vector<mType> p(d.size());
    std::vector<DivisType<cType, mType>> q(d.size());
    std::vector<DivisType<cType, mType>> x(d.size());

    p[0] = -A.c(0) / A.b(0);
    q[0] = d[0] / A.b(0);

    for(std::size_t i = 1; i < d.size() - 1; i++){
        p[i] = -(A.c(i) / (A.a(i - 1) * p[i - 1] + A.b(i)));
        q[i] = (d[i] - A.a(i - 1) * q[i - 1]) / (A.a(i - 1) * p[i - 1] + A.b(i));
    }

    q[d.size() - 1] = (d[d.size() - 1] - A.a(d.size() - 2) * q[d.size() - 2]) / (A.a(d.size() - 2) * p[d.size() - 2] + A.b(d.size() - 1));

    x[d.size() - 1] = q[d.size() - 1];

    for(std::size_t i = 1; i < d.size(); i++){
        x[d.size() - 1 - i] = p[d.size() - 1 - i] * x[d.size() - i] + q[d.size() - 1 - i];
    }

    return x;
}


/**
* xType - тип аргумента x.
* yType - тип значения функции y
*/
template<typename xType, typename yType>
class CubicSpline {
    std::vector<xType> points_;
    std::vector<yType> values_;

    using DeltaXType = DiffType<xType>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;

    struct Coeffs {
        yType a;
        yType b;
        yType c;
        yType d;
    };

    std::vector<Coeffs> coeffs_;
    std::size_t size_;
public:
    CubicSpline( const std::vector<xType> &points,  // Значения x
                 const std::vector<yType>& values,  // значения y
                 const Deriv2Type& leftDeriv,  // значение для левой второй производной
                 const Deriv2Type& rightDeriv  // значение для правой второй производной
    ):points_(points), values_(values){
        size_ = points.size();
        coeffs_.resize(size_);

        coeffs_[0].c = leftDeriv;
        coeffs_[size_ - 1].c = rightDeriv;

        std::vector<yType> u(size_ - 2);
        std::vector<xType> upThreeDiagVect(size_ - 3), centralThreeDiagVect(size_ - 2, 2), downThreeDiagVect(size_ - 3);

        std::vector<DeltaXType> h(size_);
        std::vector<DiffType<yType>> raznosti(size_);

        for(int i = 1; i < size_; i++){
            h[i] = points_[i] - points_[i - 1];
            raznosti[i] = (values_[i] - values_[i - 1]) / h[i];
        }

        u[0] = 6 * (raznosti[2] - raznosti[1]) / (points_[2] - points_[0]);
        for(int i = 0; i < size_ - 3; i++){
            upThreeDiagVect[i] = h[i + 2] / (h[i + 2] + h[i + 1]);
            downThreeDiagVect[i] = h[i + 1] / (h[i + 2] + h[i + 1]);
            u[i + 1] = 6 * (raznosti[i + 3] - raznosti[i + 2]) / (points_[i + 3] - points_[i + 1]);
        }

        ThreeDiagonalMatrix<xType> threeDiagonalMatrix(upThreeDiagVect, centralThreeDiagVect, downThreeDiagVect);
        std::vector<Deriv2Type> solution = ThreeDiagonalSolve(threeDiagonalMatrix, u);

        for(int i = 1; i < size_ - 1; i++){
            coeffs_[i].c = solution[i - 1];
            coeffs_[i].d = (coeffs_[i].c - coeffs_[i - 1].c) / h[i];
            coeffs_[i].a = values[i];
            coeffs_[i].b = h[i] / 3 * (coeffs_[i].c + coeffs_[i - 1].c / 2) + raznosti[i];
        }

        coeffs_[size_ - 1].a = values[size_ - 1];
        coeffs_[size_ - 1].b = h[size_ - 1] / 3 * (coeffs_[size_ - 1].c + coeffs_[size_ - 2].c / 2) + raznosti[size_ - 1];
        coeffs_[size_ - 1].d = (coeffs_[size_ - 1].c - coeffs_[size_ - 2].c) / h[size_ - 1];

    }

    CubicSpline( const std::vector<xType> &points,  // Значения x
                 const std::vector<yType>& values  // значения y
    ):points_(points), values_(values){
        size_ = points.size();
        coeffs_.resize(size_);

        coeffs_[0].c = 0;
        coeffs_[size_ - 1].c = 0;

        std::vector<yType> u(size_ - 2);
        std::vector<xType> upThreeDiagVect(size_ - 3), centralThreeDiagVect(size_ - 2, 2), downThreeDiagVect(size_ - 3);

        std::vector<DeltaXType> h(size_);
        std::vector<DiffType<yType>> raznosti(size_);

        for(int i = 1; i < size_; i++){
            h[i] = points_[i] - points_[i - 1];
            raznosti[i] = (values_[i] - values_[i - 1]) / h[i];
        }

        u[0] = 6 * (raznosti[2] - raznosti[1]) / (points_[2] - points_[0]);
        for(int i = 0; i < size_ - 3; i++){
            upThreeDiagVect[i] = h[i + 2] / (h[i + 2] + h[i + 1]);
            downThreeDiagVect[i] = h[i + 1] / (h[i + 2] + h[i + 1]);
            u[i + 1] = 6 * (raznosti[i + 3] - raznosti[i + 2]) / (points_[i + 3] - points_[i + 1]);
        }

        ThreeDiagonalMatrix<xType> threeDiagonalMatrix(upThreeDiagVect, centralThreeDiagVect, downThreeDiagVect);
        std::vector<Deriv2Type> solution = ThreeDiagonalSolve(threeDiagonalMatrix, u);

        for(int i = 1; i < size_ - 1; i++){
            coeffs_[i].c = solution[i - 1];
            coeffs_[i].d = (coeffs_[i].c - coeffs_[i - 1].c) / h[i];
            coeffs_[i].a = values[i];
            coeffs_[i].b = h[i] / 3 * (coeffs_[i].c + coeffs_[i - 1].c / 2) + raznosti[i];
        }

        coeffs_[size_ - 1].a = values[size_ - 1];
        coeffs_[size_ - 1].b = h[size_ - 1] / 3 * (coeffs_[size_ - 1].c + coeffs_[size_ - 2].c / 2) + raznosti[size_ - 1];
        coeffs_[size_ - 1].d = (coeffs_[size_ - 1].c - coeffs_[size_ - 2].c) / h[size_ - 1];
    }



    yType interpolate(const xType& x) const noexcept{

        std::size_t number;
        for(int i = 1; i < size_; i++){
            if(points_[i - 1] <= x && x <= points_[i]){
                number = i;
                break;
            }
        }
        DiffType<xType> dx = (x - points_[number]);
        return coeffs_[number].a + coeffs_[number].b * dx + coeffs_[number].c / 2 * dx * dx + coeffs_[number].d / 6 * dx * dx * dx;
    }
};


#endif //VYCHMATY_SPLINE_HPP
