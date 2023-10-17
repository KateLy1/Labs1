#ifndef DERIVATIVE_H_INCLUDED
#define DERIVATIVE_H_INCLUDED
#include <iostream>
#include <fstream>
#include<array>
#include <vector>
#include<cmath>

using namespace std;

// ������. ������������� � "std::cout" ��������� � ���� ������ (��� �������)
#define PrintExpression(Expression)\
     std::cout << "{" #Expression "}: " << (Expression) <<\
     		std::endl;

// ������ ������� ��� �������
template<typename Type>
void printVector(std::vector<Type>& myVector) {  // ��� �������
    for (auto element : myVector)
        std::cout << element << ' ';
    std::cout << '\n';
}

// ������ ������� ��� �������
template<typename Type>
void print2dVector(std::vector<std::vector<Type>> const &matrix) {
    for (std::vector<Type> row: matrix) {
        for (Type value: row) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}

// ������� ���������� ����������
unsigned int factorial(unsigned int n) {
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}

template<typename RealType>
struct Matrix {
    int n; // ���������� �����
    int m; // ���������� ��������
    std::vector<std::vector<RealType>> data; // ������ �������
};

// ������� ������� �������� ��������� ������� ������
template<typename RealType>
std::vector<RealType> gauss(Matrix<RealType>& A, std::vector<RealType>& b) {
// ������ ���
    for (int i = 0; i < A.n - 1; i++) {
        // ����� ������������� �������� � �������
        int maxi = i;
        RealType max_value = std::abs(A.data[i][i]);
        for (int j = i + 1; j < A.n; j++) {
            if (std::abs(A.data[j][i]) > max_value) {
                maxi = j;
                max_value = std::abs(A.data[j][i]);
            }
        }

// ������������ �����, ����� ������������ ������� ��� �� ���������
        if (maxi != i) {
            std::swap(A.data[i], A.data[maxi]);
            std::swap(b[i], b[maxi]);
        }

// �������������� ������� ������ � ����� ����
        for (int j = i + 1; j < A.n; j++) {
            RealType factor = A.data[j][i] / A.data[i][i];
            for (int k = i; k < A.m; k++) {
                A.data[j][k] -= factor * A.data[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // �������� ���
    std::vector<RealType> x(A.n);
    for (int i = A.n - 1; i >= 0; i--) {
        RealType sum = 0;
        for (int j = i + 1; j < A.m; j++) {
            sum += A.data[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A.data[i][i];
    }

    return x;
}

// ������������� ������� ��� ���������� ������������� ��� ����� points
template<typename RealType>
void InitMatrix(Matrix<RealType>& A, const std::vector<RealType>& points) {
    int n = points.size();

    for (int i = 0; i <= n; i++){       // ������� ������� (n+1)x(n+1)
        std::vector<RealType> v;
        for (int j = 0; j <= n; j++) {
            v.push_back(1);             // ������ ������ ������� ������� �� 1
        }
        A.data.push_back(v);
    }

 //   for (int j = 0; j <= n; j++)
 //       A.data[0][j] = 1;
    for (int i = 1; i <= n; i++)
        A.data[i][0] = 0;

    for (int j = 1; j <= n; j++) {
        for (int i = 1; i <= n; i++) {
            A.data[i][j] = pow(points[j-1], i)/factorial(i);
        }
    }
}

// ������������ ��� ���������� ������ �����������
template<typename RealType>
struct DerivativeCoef {
    RealType centralCoef;
    std::vector<RealType> otherCoefs;
};


// ������� ��� ���������� �������������
template<typename RealType>
DerivativeCoef<RealType> calcDerivativeCoef(const std::vector<RealType>& points) noexcept {
    DerivativeCoef<RealType> dcoef;
    Matrix<RealType> A;
    InitMatrix(A, points);
    // print2dVector<RealType>(A.data);

    int n = points.size();
    std::vector<RealType> b(n+1, 0); // ������ ����� ������� ��������� ����, ����� 2-�� ��������
    b[1] = 1;

    A.n = n + 1; A.m = n + 1;

    std::vector<RealType> x = gauss(A, b);

    dcoef.centralCoef = x[0];
    std::vector<RealType> xother(x.begin() + 1, x.end());
    dcoef.otherCoefs = xother;

    return dcoef;
}

// ���������� ������ ����������� �� �������������.
// ������������� ������� ��� ���������� � ������ (x0+h, x0-h � �.�) ���������� � ���� ���������.
template<typename RealType>
RealType calcDerivative1st(RealType func(RealType), 
                            DerivativeCoef<RealType> dcoef, 
                            const std::vector<RealType>& points, 
                            const RealType x0, 
                            const RealType h) noexcept {
    int n = points.size();
    RealType d = 0;
    for (int i = 0; i < n; i++) {
        d += dcoef.otherCoefs[i] * func(x0 + h * points[i]);
    }
    d += dcoef.centralCoef * func(x0);
    d /= h;
    return d;
}

//������ �������������
template<typename RealType>
void printCoefs(DerivativeCoef<RealType> dcoef, const std::vector<RealType>& points) {

    std::cout << "points = {";
    for (RealType p : points) {
        std::cout << p << " ";
    }
    std::cout << "}" << std::endl;

//    std::cout << points.size() << " points:" << std::endl;
    std::cout << "centralCoef = " << dcoef.centralCoef << std::endl;
    std::cout << "otherCoef = {";

    for (RealType value : dcoef.otherCoefs) {
        std::cout << value << " ";
    }

    std::cout << "}" << std::endl;
    std::cout << std::endl;
}


#endif // DERIVATIVE_H_INCLUDED
