#include "spline.hpp"

int main() {
	double a = 0, b = 10;                               // ������� [a;b]
    std::vector<int> N_nodes {5, 10, 20, 40, 80, 160};  // ��������� ���������� ����� �������
    std::vector<double> x, y;                           // x, y -���������� �����
	int np = 1000;                                      // ���������� ����� ������������
	std::vector<double> xp, yp, y_func;                 // ���������� ����� ������������ � �������� ����������� ������� � xp
	std::vector<double> splineErr;                      // ������ ��� ������ ������ ������������ ��� ������ ���������� �����
    double s = (b-a)/(np-1);
    for (int i = 0; i < np; i++) {
        xp.push_back(a + s*i);
        y_func.push_back(std::exp(a + s*i));            // �������� ������� � ������ ������������
    }

    // ���� ���������� ������� �� ������� ���������� �����
    for (auto n : N_nodes){
        double h = (b-a)/(n-1);   // ��������� �������
        for (int i = 0; i < n; i++){
            x.push_back(a + h*i);
            y.push_back(std::exp(a + h*i));
        }
        CubicSpline<double, double> spline(x, y);
        yp = spline.interpolate(xp);
        splineErr.push_back(ErrVectors<double>(yp, y_func));  // ���������� ������ ������������ � ������;

        std::vector<double>().swap(x);                        // ������� ��� �������� ������� x � ����������� ������ ��� ������������� � ������ ����������� �����
        std::vector<double>().swap(y);                        // ������� ��� �������� ������� y � ����������� ������ ��� ������������� � ������ ����������� �����
    }

    // ������ ���������� � ����
    std::ofstream output;
    output.open("splineErr.csv");
    for (unsigned int i = 0; i < N_nodes.size(); i++) {
        output << N_nodes[i] << ";" << splineErr[i] << std::endl;
    }
    output.close();

    return 0;

}

