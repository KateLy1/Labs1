#include "pch.h"
#include "../Integrate/integrate.hpp"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

double f3(double x) {
    return x * std::cos(x);
}

double f4(double x) {
    return pow(x,2) * std::cos(x);
}

double f5(double x) {
    return pow(x, 3) * std::cos(x);
}


TEST(TestCaseIntegral, TestGauss3) {

    const size_t N3 = 3;
    const double a = 0., b = M_PI; // отрезок интегрирования [a;b]
    double h = 0.1, integral;
    double integral_true = -2.;

    integral = integrate_dx<double, double, N3>(f3, a, b, h); // функцию передаем в виде параметра
    double integral_test = round(integral * 10000) / 10000;

    EXPECT_EQ(integral_true, integral_test);
}

TEST(TestCaseIntegral, TestGauss4) {

    const size_t N4 = 4;
    const double a = 0., b = M_PI; // отрезок интегрирования [a;b]
    double h = 0.1, integral;
    double integral_true = round( - 2. * M_PI * 10000) / 10000;  //  -2*PI

    integral = integrate_dx<double, double, N4>(f4, a, b, h); // функцию передаем в виде параметра
    double integral_test = round(integral * 10000) / 10000;

    EXPECT_EQ(integral_true, integral_test);
}

TEST(TestCaseIntegral, TestGauss5) {

    const size_t N5 = 5;
    const double a = 0., b = M_PI; // отрезок интегрирования [a;b]
    double h = 0.1, integral;
    double integral_true = round( - 3. * (M_PI * M_PI - 4.) * 10000) / 10000;   // - 3 * (PI^2 - 4)

    integral = integrate_dx<double, double, N5>(f5, a, b, h); // функцию передаем в виде параметра
    double integral_test = round(integral * 10000) / 10000;

    EXPECT_EQ(integral_true, integral_test);
}