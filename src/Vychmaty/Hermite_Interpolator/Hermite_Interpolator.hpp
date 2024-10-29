#ifndef VYCHMATY_HERMITE_INTERPOLATOR_HPP
#define VYCHMATY_HERMITE_INTERPOLATOR_HPP

#include <array>


template<typename xType, typename yType, unsigned int N>
class HermiteInterpolator {
    std::array<xType, 2 * N> z_;
    std::array<yType, 2 * N> f_;
    std::array<yType, N> deriv_;
    std::array<yType, 2 * N> res_Coefs_;
public:
    HermiteInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values, const std::array<yType, N>& deriv) noexcept{
        for(int i = 0; i < N; i++){
            z_[2 * i] = points[i];
            z_[2 * i + 1] = points[i];

            f_[2 * i] = values[i];
            f_[2 * i + 1] = values[i];

            deriv_[i] = deriv[i];
        }

        res_Coefs_ = f_;
        for(int i = 1; i < 2 * N; i++){
            for(int j = 2 * N - 1; j >= i; j--){
                if(z_[j] == z_[j - i]){
                    res_Coefs_[j] = deriv_[(j - i) / 2];
                } else{
                    res_Coefs_[j] = (res_Coefs_[j] - res_Coefs_[j - 1]) / (z_[j] - z_[j - i]);
                }
            }
        }
    }

    yType interpolate(const xType& x) const noexcept{
        yType res, tmp;

        res = res_Coefs_[0];
        for(int i = 1; i < 2 * N; i++){
            tmp = res_Coefs_[i];
            for(int j = 0; j < i; j++){
                tmp *= (x - z_[j]);
            }
            res += tmp;
        }

        return res;
    }
};


#endif //VYCHMATY_HERMITE_INTERPOLATOR_HPP
