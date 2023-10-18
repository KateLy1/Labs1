#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <fstream>

using namespace std;

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

template <typename RealType>
struct GaussQ {
    const std::vector<RealType> xg2{-1./sqrt(3.), 1./sqrt(3.)};
    const std::vector<RealType> wg2{1., 1.};
    const std::vector<RealType> xg3{-sqrt(3./5.), 0., sqrt(3./5.)};
    const std::vector<RealType> wg3{5./9., 8./9., 5./9.};
    const std::vector<RealType> xg4{-sqrt(3./7. + 2./7.*sqrt(6./5.)), -sqrt(3./7.-2./7.*sqrt(6./5.)), sqrt(3./7.-2./7.*sqrt(6./5.)), sqrt(3./7. + 2./7.*sqrt(6./5.))};
    const std::vector<RealType> wg4{(18.-sqrt(30.))/36., (18.+sqrt(30.))/36., (18.+sqrt(30.))/36., (18.-sqrt(30.))/36.};
    const std::vector<RealType> xg5{-1./3.*sqrt(5.+2.*sqrt(10./7.)),-1./3.*sqrt(5.-2.*sqrt(10./7.)), 0., 1./3.*sqrt(5.-2.*sqrt(10./7.)), 1./3.*sqrt(5.+2.*sqrt(10./7.))};
    const std::vector<RealType> wg5{(322.-13.*sqrt(70.))/900., (322.+13.*sqrt(70.))/900., 128./225., (322.+13.*sqrt(70.))/900., (322.-13.*sqrt(70.))/900.};
};


// Функция интегрирования на одном отрезке
template <typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate(
                         Callable func(RealType),
                         const RealType& a,
                         const RealType& b) {
    const GaussQ<RealType> polyn;
    std::vector<RealType> xg;
    std::vector<RealType> wg;

    switch ( N )  // выбор квадратуры Гаусса N от 2 до 5
      {
         case 2:
            xg = polyn.xg2;
            wg = polyn.wg2;
            break;
         case 3:
            xg = polyn.xg3;
            wg = polyn.wg3;
            break;
         case 4:
            xg = polyn.xg4;
            wg = polyn.wg4;
            break;
         case 5:
            xg = polyn.xg5;
            wg = polyn.wg5;
            break;
      }

    RealType bma = (b - a) / 2;
    RealType bpa = (b + a) / 2;

    RealType sum = 0;
    for (std::size_t i = 0; i < N; i++) {
        sum += wg[i] * func(bma * xg[i] + bpa);
    }

    return bma * sum;
}


// Функция разбивает отрезок [a;b] на подотрезки не более dx и призводит интегрирование на всех подотрезках
template <typename Callable, typename RealType, std::size_t N>
decltype(auto) integrate_dx(
                         Callable func(RealType),
                         const RealType& a,
                         const RealType& b,
                         const Dif<RealType>& dx) {

    const GaussQ<RealType> polyn;

    RealType integral = 0, integral_dx = 0, an, bn;
    int ndx = static_cast<int>((b - a)/dx);         // количество целых интервалов dx

    for (int i = 0; i <= ndx; i++) {

        if (i == ndx)
            bn = b;
        else
            bn = a + (i+1)*dx;

        an = a + i*dx;

        integral = integrate<Callable, RealType, N>(func, an, bn);

        integral_dx += integral;
    }

    return integral_dx;
}

#endif // INTEGRATE_H_INCLUDED
