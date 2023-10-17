#include "newton_interpolator.hpp"
using namespace std;

int main()
{
    // ����������� ������ ������ cout � ���� csv
    std::ofstream out("newton_interpolation.csv");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf());               //redirect std::cout to file
    std::ofstream output;

//    const std::array<int, 3> N_nodes{3, 4, 5};
    double a = 0.0, b = 2.0, bt, err;          // ������� ������������ [a;b]


// ���������� ������ ������������ ��� 3, 4 ��� 5 ����� ������������
// �� ������������������ �������� [0,2], [0,1], [0,1/2], ... [0, 1/16].
// ����������� ������������ �����.
    bool Cheb_nodes = false;
    bt = b;
    for (int bk = 0; bk<6; bk++){           // ���� �� �������� [a;b] (b �������� �� 2 �� 1/16)
        err = Newton_interpolator_err<double, double, 3>(a, bt, Cheb_nodes);
        std::cout << "3" << ";" << bt << ";" << err << ";" << "����" << std::endl;
        bt = bt/2.0;
    }
    bt = b;
    for (int bk = 0; bk<6; bk++){
        err = Newton_interpolator_err<double, double, 4>(a, bt, Cheb_nodes);
        std::cout << "4" << ";" << bt << ";" << err << ";" << "����" << std::endl;
        bt = bt/2.0;
    }
    bt = b;
    for (int bk = 0; bk<6; bk++){
        err = Newton_interpolator_err<double, double, 5>(a, bt, Cheb_nodes);
        std::cout << "5" << ";" << bt << ";" << err << ";" << "����" << std::endl;
        bt = bt/2.0;
    }

// ���������� ������ ������������ ��� 3, 4 ��� 5 ����� ������������
// �� ������������������ �������� [0,2], [0,1], [0,1/2], ... [0, 1/16].
// ����������� ������������ �����.
    Cheb_nodes = true;
    bt = b;
    for (int bk = 0; bk<6; bk++){
        err = Newton_interpolator_err<double, double, 3>(a, bt, Cheb_nodes);
        std::cout << "3" << ";" << bt << ";" << err << ";" << "���" << std::endl;
        bt = bt/2.0;
    }
    bt = b;
    for (int bk = 0; bk<6; bk++){
        err = Newton_interpolator_err<double, double, 4>(a, bt, Cheb_nodes);
        std::cout << "4" << ";" << bt << ";" << err << ";" << "���" << std::endl;
        bt = bt/2.0;
    }
    bt = b;
    for (int bk = 0; bk<6; bk++){
        err = Newton_interpolator_err<double, double, 5>(a, bt, Cheb_nodes);
        std::cout << "5" << ";" << bt << ";" << err << ";" << "���" << std::endl;
        bt = bt/2.0;
    }


    std::cout.rdbuf(coutbuf); //reset to standard output again


    return 0;
}


