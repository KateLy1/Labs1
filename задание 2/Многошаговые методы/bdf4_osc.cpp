#include <iostream>
#include <array>
#include <vector>
#include<cmath>
#include <fstream>
#include <iomanip>


double func_sin(double x)
{
    return std::sin(x);
}

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

struct RK4Table{        // ������� ������� ��� �����-����� 4-�� �������
    static constexpr unsigned int stages = 4;
    std::array<std::array<double, stages>, stages> table = {{ {0., 0., 0., 0.}, {1./2., 0., 0., 0.}, {0., 1./2., 0., 0.}, {0., 0., 1., 0.} }};
    std::array<double, stages> cColumn = {0., 1./2. , 1./2., 1.};
    std::array<double, stages> bString = {1./6., 1./3., 1./3., 1./6.};
};

/* ������������ 4-� �������� ������ */
struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size> alpha = {-3./25., 16./25., -36./25, 48./25.};
    static constexpr double beta = 12./25.;
};

struct IntegrationParameters{
    double step;  // ��� ��������������
    double epsilon; // �������� ������� ����������� ���������
    double maxIter; // ������������ ���������� �������� ��� ������� ����������� ���������
};

/*  ����� ������ ����� ����.��������� y'' = -y */
class Oscillator {
public:
        //static constexpr unsigned int dim = 2;
        using Argument = double;  //��� x
        using State = std::vector<double>;

        struct StateAndArg {
            State state {0., 0.};  // ������������� 2-������� �������
            Argument arg;
        } saa;

        // ��������� ������� ���� ��������� (������� ����������� ���������)
        //double calcDif(const State& first, const State& second) {}

        // ��������� ������ ����� �� - ������� f
        std::vector<double> calc(const StateAndArg& stateAndArg) const {
            return std::vector<double>{stateAndArg.state[1], -stateAndArg.state[0]};
        }
};

template<typename Table, typename RHS>
std::vector<double> integrate(
                      const typename RHS::StateAndArg& initialState,
                      IntegrationParameters& parameters,
                      unsigned int n)
{
    Table tableB;
    BDF4 bd;
    double A4;  // ��������������� ����������
    std::vector<double> x(n, 0), y(n, 0.), v(n, 0);  // ������������� ������
    RHS rhs;
    double step = parameters.step;
    unsigned int srk = sizeof(tableB.table[0])/sizeof(double);   // ������� ������ �����-�����
    unsigned int sbdf = bd.size;                                 // ������� ��� ������

    x[0] = initialState.arg;
    y[0] = initialState.state[0];
    v[0] = initialState.state[1];

    unsigned int s = sizeof(tableB.table[0])/sizeof(double);   // ����������� ������ �����-�����
    std::vector<double>  Ks(srk, 0), Ps(s, 0);
    std::vector<double> yv(2, 0);

    // ������ (����������� ��������� �������� y[i]) ������� �����-�����
    for (unsigned int m = 1; m < sbdf; m++) {   // ���� �� ������ 1-3 ��� �������
        x[m] = x[0] + m*step;
        y[m] = y[m-1];
        v[m] = v[m-1];
        rhs.saa.state[0] = y[m-1];
        rhs.saa.state[1] = v[m-1];
        yv = rhs.calc(rhs.saa);
        Ks[0] = yv[0];
        Ps[0] = yv[1];

        for(unsigned int j = 1; j < srk; j++) {
            rhs.saa.state[0] = y[m-1] + Ks[j-1]*tableB.table[j][j-1]*step;
            rhs.saa.state[1] = v[m-1] + Ps[j-1]*tableB.table[j][j-1]*step;
            yv = rhs.calc(rhs.saa);
            Ks[j] = yv[0];
            Ps[j] = yv[1];
        }

        for(unsigned int j = 0; j < sbdf; j++) {
            y[m] = y[m] + Ks[j]*tableB.bString[j]*step;
            v[m] = v[m] + Ps[j]*tableB.bString[j]*step;
        }
    }
    // ����� �������

    // ������������ ����� ���
    for (unsigned int m = sbdf; m < n; m++) {   // ���� �� ������ ��������� �������
        x[m] = x[0] + m*step;
        rhs.saa.arg = x[m];
        A4 = bd.alpha[0] * (y[m-4] + bd.beta*step*v[m-4])
           + bd.alpha[1] * (y[m-3] + bd.beta*step*v[m-3])
           + bd.alpha[2] * (y[m-2] + bd.beta*step*v[m-2])
           + bd.alpha[3] * (y[m-1] + bd.beta*step*v[m-1])
        ;
        y[m] = A4 / (1. + std::pow(bd.beta*step, 2));

        v[m] = bd.alpha[0] * v[m-4]
             + bd.alpha[1] * v[m-3]
             + bd.alpha[2] * v[m-2]
             + bd.alpha[3] * v[m-1]
             - bd.beta*step*y[m];
    }

    // ����� ������� ��� ��������
    /*
    std::ofstream output;
    output.open("BDF4_st.csv");
    for (unsigned int i = 0; i < n; i++){
        output << x[i] << ";" << y[i] << ";" << v[i] << std::endl;
    }
    */

    return y;
}

int main()
{
    double a = 0., b = 5., yInitialOsc = 0., vInitialOsc = 1.;
    std::vector<double> steps {1, 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6};    // ���� ��������������
    IntegrationParameters params;

    Oscillator oInitial;
    oInitial.saa.arg = a;
    oInitial.saa.state[0] = yInitialOsc;
    oInitial.saa.state[1] = vInitialOsc;

    // ������ ������ � ����
    std::ofstream output;
    output.open("BDF4oscErr.csv");
    output << std::setprecision(10);


    // ���� ���������� ������� ��� ������ �����
    for (auto step : steps){
        params.step = step;
        unsigned int m = int((b-a)/params.step) + 1;                // ���������� ����� �������
        std::vector<double> x(m, 0), y(m, 0);           // x,y ������� �������
        std::vector<double> y_solution(m, 0);           // ������ ������������ �������, ��������������� ������
        for (unsigned int i = 0; i < m; i++){
            x[i] = (a + step*i);
            y[i] = func_sin(x[i]);
        }
        y_solution = integrate<RK4Table, Oscillator>(oInitial.saa, params, m);    // ���������� ������� �������
        output << step << ";" << ErrVectors<double>(y, y_solution) << std::endl;  // ���������� ������ � ����;
    }

    output.close();
    return 0;
}
