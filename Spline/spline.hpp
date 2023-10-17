#ifndef SPLINE_HPP_INCLUDED
#define SPLINE_HPP_INCLUDED

#include<iostream>
#include<vector>
#include<cmath>
#include <fstream>

using namespace std;

// ������. ������������� � "std::cout" ��������� � ���� ������ (��� �������)
#define PrintExpression(Expression)\
     std::cout << "{" #Expression "}: " << (Expression) <<\
     		std::endl;

template<typename Type>
void printVector(std::vector<Type>& myVector)  // ��� �������
{
    for (auto element : myVector)
        std::cout << element << ' ';
    std::cout << '\n';
}


template<typename Type>
class ThreeDiagonalMatrix
    // ���������� 3-� ���������� ������� � ������ ����� R ��� ������������� �������
    // n - ����� ��������� (����� �������)
    // b - ���������, ������� ��� ������� (����������: [0;n-2]) �.�. b[n-1] = 0
    // c - ������� ��������� ������� (����������: [0;n-1])
    // a - ���������, ������� ��� ������� (����������: [1;n-1]) �.�. a[0] = 0
    // R - ������ ����� ������� ���������   R[0] = R[N-1] = 0 - ������������ ������ (��������� �������)
{
    std::vector<Type> dx, dy;
public:
    ThreeDiagonalMatrix(const std::vector<Type>& dpoints,
        const std::vector<Type>& dvalues) : dx(dpoints), dy(dvalues)
    {
    }
    void m_init(std::vector<Type>& adiag, std::vector<Type>& cdiag, std::vector<Type>& bdiag, std::vector<Type>& R)
    {
        int n = dx.size() + 1;      // ����������� ���������� ����� ���������� ����� - �� 1 ������, ��� ���������� dx
        cdiag.push_back(1);         // c[0] = 1
        bdiag.push_back(0);         // b[0] = 0
        adiag.push_back(0);         // a[0] = 0
        R.push_back(0);             // R[0] = 0     // ������ �����
        for (int i = 1; i < n - 1; i++) {
            adiag.push_back(dx[i - 1]);
            cdiag.push_back(2 * (dx[i - 1] + dx[i]));
            bdiag.push_back(dx[i]);
            R.push_back(3 * (dy[i] / dx[i] - dy[i - 1] / dx[i - 1]));
        }
        cdiag.push_back(1);         // c[n-1] = 1;
        bdiag.push_back(0);         // b[n-1] = 0
        adiag.push_back(0);         // a[0] = 0
        R.push_back(0);             // R[n-1]
    }

    std::vector<Type> adiag()
    {
        std::vector<Type> ad, cd, bd, R;
        m_init(ad, cd, bd, R);
        return ad;
    }

    std::vector<Type> cdiag()
    {
        std::vector<Type> ad, cd, bd, R;
        m_init(ad, cd, bd, R);
        return cd;
    }

    std::vector<Type> bdiag()
    {
        std::vector<Type> ad, cd, bd, R;
        m_init(ad, cd, bd, R);
        return bd;
    }

    std::vector<Type> F()     // ������ �����
    {
        std::vector<Type> ad, cd, bd, R;
        m_init(ad, cd, bd, R);
        return R;
    }
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

// ������� ��� ������� ������� ��������. x - ������� (������)
// ������� ������ �������: a[0] = 0, b[n-1] = 0
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>>solve(ThreeDiagonalMatrix<mType>& matrix,
                                          std::vector<cType>& f)
{
    std::vector<mType> a, c, b;  // ���������� ����������.
    c = matrix.cdiag();
    b = matrix.bdiag();
    a = matrix.adiag();
    //    f = matrix.F();           // f - ������ ����� - ����� ������ �����, � ����� �������� � ���������� �������

    // ����� ��������
    mType m;
    int n = a.size();
    std::vector<mType> x(n);
    for (int i = 1; i < n; i++)
    {
        m = a[i] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    x[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
    }
    return x;

}

// ����� Spline ��������
//   ������� �� ���������� ������������� ���������� ��������� ��� ������� ��������� ��������� �������
//   � ������� ������������ ��� ������ ����� �� �������
template<typename xType, typename yType>
class CubicSpline
{
    std::vector<xType> xn;
    std::vector<yType> yn;
public:
    CubicSpline(const std::vector<xType>& points,
        const std::vector<yType>& values) : xn(points), yn(values)   // ���� �������
    {
    }
    // ������� ���������� ������������� �������� �� ������ ���������
    void SplineCoeff(std::vector<xType>& ai,
        std::vector<xType>& bi,
        std::vector<xType>& ci,
        std::vector<xType>& di,
        std::vector<xType>& M)
    {
        int n = int(xn.size());                 // ���������� ����� �������

        // ����������� �������� ������������
        std::vector<xType> dx, dy;
        for (int i = 0; i < n - 1; i++) {         // ���������� ���������� ������ �� 1 ���������� �����
            dx.push_back(xn[i + 1] - xn[i]);
            dy.push_back(yn[i + 1] - yn[i]);
        }

        // ������� ������� ���������
        std::vector<xType> F;
        ThreeDiagonalMatrix<xType> diagMat(dx, dy);     //  ������������� 3-� ������������ �������
        F = diagMat.F();                         // ������������� ������ ����� �������� � ����� 3� �������
        M = solve(diagMat, F);                   //  ������� ������� M
        //printVector(M);

    // ���������� ������������ ����������� �������� ��� ������� ���������
        for (int i = 0; i < n - 1; i++) {
            ai.push_back(yn[i]);
            di.push_back((M[i + 1] - M[i]) / (3 * dx[i]));
            bi.push_back(dy[i] / dx[i] - dx[i] * (2 * M[i] + M[i + 1]) / 3);
            ci.push_back(M[i]);
        }
    }

    // ������� ������������ ��� ������ ����� (�������)
    std::vector<yType> interpolate(const std::vector<xType>& xp)
    {
        std::vector<xType> ai, bi, ci, di, M;
        std::vector<yType> yp;
        SplineCoeff(ai, bi, ci, di, M);
        int n = int(xn.size());                 // ���������� ����� �������
        int np = int(xp.size());                // ���������� ����� ������������

        // ���������� ������� � ������ ������������
        for (int i = 0; i < np; i++) {
            int k = 0;
            for (int j = 0; j < n - 1; j++) {
                if (xp[i] >= xn[j] && xp[i] < xn[j + 1]) {
                    k = j;
                    break;
                }
                else if (xp[i] == xn[n - 1]) {
                    k = n - 2;
                }
            }
            yType temp = ai[k] + bi[k] * (xp[i] - xn[k]) + M[k] * pow((xp[i] - xn[k]), 2) + di[k] * pow((xp[i] - xn[k]), 3);
            yp.push_back(temp);
        }
        return yp;
    }
};
// ������� ���������� ������������ ������� �� ������ ��������� ���� ��������
template <typename yType>
yType ErrVectors(std::vector<yType>& v1, std::vector<yType>& v2)
{
    int n = int(v1.size());
    yType err_mod_max = 0.0, err_mod = 0.0;
    for (int i = 0; i < n; i++) {
        err_mod = std::abs(v1[i] - v2[i]);
        if (err_mod > err_mod_max)
            err_mod_max = err_mod;
    }
    return err_mod_max;
}

#endif // SPLINE_HPP_INCLUDED
