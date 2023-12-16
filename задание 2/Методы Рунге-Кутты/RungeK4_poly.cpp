#include <iostream>
#include <array>
#include <vector>
#include<cmath>
#include <fstream>
#include <iomanip>


double func(double x)  // функция аналитического решения
{
    return pow(x,4)/4.;
}


// функция вычисления максимальной разницы по модулю координат двух векторов
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

struct RK4Table{        // таблица Бутчера для Рунге-Кутты 4-го порядка
    static constexpr unsigned int stages = 4;
    std::array<std::array<double, stages>, stages> table = {{ {0., 0., 0., 0.}, {1./2., 0., 0., 0.}, {0., 1./2., 0., 0.}, {0., 0., 1., 0.} }};
    std::array<double, stages> cColumn = {0., 1./2. , 1./2., 1.};
    std::array<double, stages> bString = {1./6., 1./3., 1./3., 1./6.};
};

class Polynom {
public:
        using Argument = double;  //тип x
        using State = double;

        struct StateAndArg {
            State state;
            Argument arg;
        } saa;

        double calc(const StateAndArg& stateAndArg) const {
            return pow(stateAndArg.arg, 3);
        }
};

template<typename Table, typename RHS>
std::vector<double> integrate(
                      const typename RHS::StateAndArg& initialState,
                      double step,
                      int n)
{
    Table tableB;
    std::vector<double> x(n, 0.), y(n, 0.);  // инициализация нулями
    RHS rhs;

    x[0] = initialState.arg;
    y[0] = initialState.state;   //y0;
    unsigned int s = sizeof(tableB.table[0])/sizeof(double);   // размерность метода Рунге-Кутты
    std::vector<double>  Ks(s, 0);

    for (int m = 1; m < n; m++) {
        x[m] = x[0] + m*step;
        y[m] = y[m-1];
        rhs.saa.arg = x[m-1] + tableB.cColumn[0]*step;
        Ks[0] = rhs.calc(rhs.saa);
        for(unsigned int j = 1; j < s; j++) {
            rhs.saa.arg = x[m-1] + tableB.cColumn[j]*step;
            Ks[j] = rhs.calc(rhs.saa);
        }
        for(unsigned int j = 0; j < s; j++) {
            y[m] = y[m] + Ks[j]*tableB.bString[j]*step;
        }
    }

    //for (int i = 0; i < n; i++){
    //    std::cout << x[i] << ";" << y[i] << ";" << std::endl;
    // }

    return y;
}

int main()
{
    double a = 0., b = 5., yInitialPolynom = 0.;     // отрезок на котором находится решение и начальное условие
    std::vector<double> steps {1, 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6};    // шаги интегрирования
    Polynom pInitial;
    pInitial.saa.arg = a;
    pInitial.saa.state = yInitialPolynom;

    // запись ошибки в файл
    std::ofstream output;
    output.open("polynomErr.csv");
    output << std::setprecision(10);

    // цикл нахождения решения для разных шагов
    for (auto step : steps){
        int m = int((b-a)/step) + 1;                // количество узлов решения
        std::vector<double> x(m, 0), y(m, 0);           // x,y точного решения
        std::vector<double> y_solution(m, 0);           // вектор вычисленного решения, инициализирован нулями
        for (int i = 0; i < m; i++){
            x[i] = (a + step*i);
            y[i] = func(x[i]);
        }
        y_solution = integrate<RK4Table, Polynom>(pInitial.saa, step, m);         // нахождение вектора решения
        output << step << ";" << ErrVectors<double>(y, y_solution) << std::endl;  // записываем ошибку в файл;
    }

    output.close();

    return 0;
}
