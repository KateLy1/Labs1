#include "pch.h"
#include "../Spline/spline.hpp"

TEST(TestCaseSpline, TestLinear) {    // тест интерпол€ции сплайном линейной функции y = x в точке (xtest, ytest)
  double a = 0, b = 5;                // отрезок [a;b]
  double xtest = 3.7, ytest, yfunc;   
  int n = 101;                        // количество узлов сплайна дл€ теста
  std::vector<double> x, y;           // x, y -координаты узлов
  std::vector<double> xp, yp;         // тестова€ точка интерпол€ции, заданна€ вектором (в тесте размерность 1)
  xp.push_back(xtest);

  double h = (b - a) / (n - 1);   
  for (int i = 0; i < n; i++) {   
       x.push_back(a + h * i);
       y.push_back(a + h * i);
  }
  CubicSpline<double, double> spline(x, y);
  yp = spline.interpolate(xp);                      // вычисление интерпол€ционного значени€ в точке xtest
  ytest = round(yp[0] * 1000) / 1000;              // округление интерпол€ции до 3х знаков после зап€той
  yfunc = round(xtest * 1000) / 1000;              // вычисление значени€ функции в точке xtest с округлением

  EXPECT_EQ(yfunc, ytest);                         //значение функции должно совпадать с интерпол€цией с точностью до округлени€
}

TEST(TestCaseSpline, TestCosinus) {     // тест интерпол€ции сплайном функции y = сos(x) в точке (xtest, ytest)
    double a = 0, b = 5;                // отрезок [a;b]
    double xtest = 3.7, ytest, yfunc;
    int n = 11;                        // количество узлов сплайна дл€ теста
    std::vector<double> x, y;          // x, y -координаты узлов
    std::vector<double> xp, yp;        // тестова€ точка интерпол€ции, заданна€ вектором (в тесте размерность 1)
    xp.push_back(xtest);

    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        x.push_back(a + h * i);
        y.push_back(std::cos(a + h * i));
    }
    CubicSpline<double, double> spline(x, y);
    yp = spline.interpolate(xp);                      // вычисление интерпол€ционного значени€ в точке xtest
    ytest = round(yp[0] * 1000) / 1000;               // округление до 3х знаков после зап€той
    yfunc = round(std::cos(xtest) * 1000) / 1000;     // вычисление значени€ функции в точке xtest с округлением

    EXPECT_EQ(yfunc, ytest);                          //значение функции должно совпадать с интерпол€цией с точностью до округлени€
}