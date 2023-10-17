#ifndef NEWTON_INTERPOLATOR_HPP_INCLUDED
#define NEWTON_INTERPOLATOR_HPP_INCLUDED
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// ����� ������������ �������
template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator
{
std::array<xType, N> Xm; // ���� ������������
std::array<yType, N> Ym;

public:
    NewtonInterpolator(const std::array<xType, N>& points,
                       const std::array<yType, N>& values) noexcept : Xm(points), Ym(values)
    {
    }

    yType interpolate(const xType& x) const noexcept   // x - ����� ������������
    {
        yType DN, f, denom;  // DN-�������, denom-����������� ����������� ��������
        unsigned int i,j,k;
        DN = Ym[0];
        for(i=1; i<N; i++){
            f = 0;
            //��������� ��������� ��������
            for(j=0; j<=i; j++){
                denom = 1;
                //������� ����������� ����������� ��������
                for(k=0; k<=i; k++){
                    if (k!=j) denom *= (Xm[j]-Xm[k]);
                }
                //������� ����������� ��������
                f += Ym[j]/denom;
            }
            //��������� ����������� �������� �� ������ (x-x[0])...(x-x[i-1])
            for(k=0; k<i; k++)
                f *= (x-Xm[k]);
            DN += f;   //�������
        }
    return DN;
    }
};

// ������� ������������ ������� ������� �� ���� ������������ (����� ����� N>1)
template <typename xType, unsigned int N>
std::array<xType, N> divide_ab(xType a, xType b, std::array<xType, N>& x_nodes)
{
    int nn = int(x_nodes.size());
    xType h = (b-a)/(nn-1);
    for (int i = 0; i < nn; i++) {
        x_nodes[i] = a + h*i;
    }
    return x_nodes;
}

// ������� ������� ������� �� ����������� ���� ������������ (����� ����� N>1)
template <typename xType, unsigned int N>
std::array<xType, N> divide_ab_Cheb(xType a, xType b, std::array<xType, N>& x_nodes)
{
    int nn = int(x_nodes.size());
    for (int i = 0; i < nn; i++) {
        x_nodes[nn-1-i] = (b+a)/2 + ((b-a)/2) * std::cos(M_PI*(2*i+1)/(2*nn));  // ���������� ������� ������, ����� ���� ��� � ������������ �������
    }
    return x_nodes;
}


// ������� ���������� ������ ������������. N_points - ���������� ����� ������������
template<typename xType, typename yType, unsigned int N_points>
yType Newton_interpolator_err(xType a, xType b, bool Cheb_nodes)    // ����������� ��� ����������� ����
{
    yType y_interpolated, y_true;
    std::array<xType, N_points> Xm {};
    std::array<yType, N_points> Ym {};
    yType err_mod_max, err_mod;
    const unsigned int ni_points = 1000;    // ���������� ����� ������������
    std::array<xType, ni_points> Xi {};    // ����� ������������

    err_mod_max = 0.0; err_mod = 0.0;
    std::fill(std::begin(Xm), std::end(Xm), 0.0);  // ��������� ��������
    std::fill(std::begin(Ym), std::end(Ym), 0.0);
    std::fill(std::begin(Xi), std::end(Xi), 0.0);

    if (Cheb_nodes)
        divide_ab_Cheb<xType, N_points> (a, b, Xm);   // ������� ������� �� ����� ������������
    else
        divide_ab<xType, N_points> (a, b, Xm);

    for (unsigned int i=0; i<N_points; i++){  // ���������� ������� � �����
        Ym[i] = std::exp(Xm[i]);
    }
    NewtonInterpolator<xType, yType, N_points> Ni(Xm, Ym);  // ���������� ���������� �� �����

    divide_ab<xType, ni_points>(a, b, Xi);     //������� ������� �� ������ ������������
    for (unsigned int i=0; i<ni_points; i++) {  // ���� �� ���� ������ ������������ � ���������� ������������ ������
        y_interpolated = Ni.interpolate(Xi[i]);      // ���������� ����������������� �������� ������� � �����
        y_true = std::exp(Xi[i]);                   // ���������� ������� �������� ������� �� �������
        err_mod = std::abs(y_true - y_interpolated);
        if (err_mod > err_mod_max) err_mod_max = err_mod;
    }

    return err_mod_max;
}

#endif // NEWTON_INTERPOLATOR_HPP_INCLUDED
