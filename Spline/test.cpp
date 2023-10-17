#include "pch.h"
#include "../Spline/spline.hpp"

TEST(TestCaseSpline, TestLinear) {    // ���� ������������ �������� �������� ������� y = x � ����� (xtest, ytest)
  double a = 0, b = 5;                // ������� [a;b]
  double xtest = 3.7, ytest, yfunc;   
  int n = 101;                        // ���������� ����� ������� ��� �����
  std::vector<double> x, y;           // x, y -���������� �����
  std::vector<double> xp, yp;         // �������� ����� ������������, �������� �������� (� ����� ����������� 1)
  xp.push_back(xtest);

  double h = (b - a) / (n - 1);   
  for (int i = 0; i < n; i++) {   
       x.push_back(a + h * i);
       y.push_back(a + h * i);
  }
  CubicSpline<double, double> spline(x, y);
  yp = spline.interpolate(xp);                      // ���������� ����������������� �������� � ����� xtest
  ytest = round(yp[0] * 1000) / 1000;              // ���������� ������������ �� 3� ������ ����� �������
  yfunc = round(xtest * 1000) / 1000;              // ���������� �������� ������� � ����� xtest � �����������

  EXPECT_EQ(yfunc, ytest);                         //�������� ������� ������ ��������� � ������������� � ��������� �� ����������
}

TEST(TestCaseSpline, TestCosinus) {     // ���� ������������ �������� ������� y = �os(x) � ����� (xtest, ytest)
    double a = 0, b = 5;                // ������� [a;b]
    double xtest = 3.7, ytest, yfunc;
    int n = 11;                        // ���������� ����� ������� ��� �����
    std::vector<double> x, y;          // x, y -���������� �����
    std::vector<double> xp, yp;        // �������� ����� ������������, �������� �������� (� ����� ����������� 1)
    xp.push_back(xtest);

    double h = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) {
        x.push_back(a + h * i);
        y.push_back(std::cos(a + h * i));
    }
    CubicSpline<double, double> spline(x, y);
    yp = spline.interpolate(xp);                      // ���������� ����������������� �������� � ����� xtest
    ytest = round(yp[0] * 1000) / 1000;               // ���������� �� 3� ������ ����� �������
    yfunc = round(std::cos(xtest) * 1000) / 1000;     // ���������� �������� ������� � ����� xtest � �����������

    EXPECT_EQ(yfunc, ytest);                          //�������� ������� ������ ��������� � ������������� � ��������� �� ����������
}