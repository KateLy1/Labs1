#include "pch.h"
#include "../Interpol/newton_interpolator.hpp"

TEST(Newton_interpolation, Uniform) {             // тест интерпол€ции функции y = x^3. –авномерные узлы интерпол€ции.
  const unsigned int N_points = 6;
  double a = 1, b = 3;
  double xtest = 1.9, ytest, yfunc;
  std::array<double, N_points> Xm{};
  std::array<double, N_points> Ym{};

  divide_ab<double, N_points>(a, b, Xm);                    
  for (unsigned int i = 0; i < N_points; i++) {             
      Ym[i] = pow(Xm[i], 3);
  }

  NewtonInterpolator<double, double, N_points> Ni(Xm, Ym);  // построение многочлена по узлам
  ytest = Ni.interpolate(xtest);                            // вычисление интерпол€ционного значени€ функции в точке xtest
  ytest = round(ytest * 1000) / 1000;                      // округление до 3х знаков после зап€той
  yfunc = round(pow(xtest, 3) * 1000) / 1000;              // вычисление значени€ функции в точке xtest

  EXPECT_EQ(yfunc, ytest);
}

TEST(Newton_interpolation, Chebyshev) {             // тест интерпол€ции функции y = x^3. „ебышевские узлы интерпол€ции.
    const unsigned int N_points = 6;
    double a = 1, b = 3;
    double xtest = 1.9, ytest, yfunc;
    std::array<double, N_points> Xm{};
    std::array<double, N_points> Ym{};

    divide_ab_Cheb<double, N_points>(a, b, Xm);                    
    for (unsigned int i = 0; i < N_points; i++) {             
        Ym[i] = pow(Xm[i], 3);
    }

    NewtonInterpolator<double, double, N_points> Ni(Xm, Ym);  // построение многочлена по узлам
    ytest = Ni.interpolate(xtest);                            // вычисление интерпол€ционного значени€ функции в точке xtest
    ytest = round(ytest * 1000) / 1000;                      // округление до 3х знаков после зап€той
    yfunc = round(pow(xtest, 3) * 1000) / 1000;              // вычисление значени€ функции в точке xtest

    EXPECT_EQ(yfunc, ytest);
}