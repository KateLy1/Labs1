#include "pch.h"
#include "../Interpol/newton_interpolator.hpp"

TEST(Newton_interpolation, Uniform) {             // ���� ������������ ������� y = x^3. ����������� ���� ������������.
  const unsigned int N_points = 6;
  double a = 1, b = 3;
  double xtest = 1.9, ytest, yfunc;
  std::array<double, N_points> Xm{};
  std::array<double, N_points> Ym{};

  divide_ab<double, N_points>(a, b, Xm);                    
  for (unsigned int i = 0; i < N_points; i++) {             
      Ym[i] = pow(Xm[i], 3);
  }

  NewtonInterpolator<double, double, N_points> Ni(Xm, Ym);  // ���������� ���������� �� �����
  ytest = Ni.interpolate(xtest);                            // ���������� ����������������� �������� ������� � ����� xtest
  ytest = round(ytest * 1000) / 1000;                      // ���������� �� 3� ������ ����� �������
  yfunc = round(pow(xtest, 3) * 1000) / 1000;              // ���������� �������� ������� � ����� xtest

  EXPECT_EQ(yfunc, ytest);
}

TEST(Newton_interpolation, Chebyshev) {             // ���� ������������ ������� y = x^3. ����������� ���� ������������.
    const unsigned int N_points = 6;
    double a = 1, b = 3;
    double xtest = 1.9, ytest, yfunc;
    std::array<double, N_points> Xm{};
    std::array<double, N_points> Ym{};

    divide_ab_Cheb<double, N_points>(a, b, Xm);                    
    for (unsigned int i = 0; i < N_points; i++) {             
        Ym[i] = pow(Xm[i], 3);
    }

    NewtonInterpolator<double, double, N_points> Ni(Xm, Ym);  // ���������� ���������� �� �����
    ytest = Ni.interpolate(xtest);                            // ���������� ����������������� �������� ������� � ����� xtest
    ytest = round(ytest * 1000) / 1000;                      // ���������� �� 3� ������ ����� �������
    yfunc = round(pow(xtest, 3) * 1000) / 1000;              // ���������� �������� ������� � ����� xtest

    EXPECT_EQ(yfunc, ytest);
}