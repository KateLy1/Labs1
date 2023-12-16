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

/* Коэффициенты 4-х шагового метода */
struct BDF4{
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size> alpha = {-3./25., 16./25., -36./25, 48./25.};
    static constexpr double beta = 12./25.;
};

struct IntegrationParameters{
    double step;  // шаг интегрирования
    double epsilon; // точность решения нелинейного уравнения
    double maxIter; // максимальное количество итераций для решения нелинейного уравнения
};

/*  Класс правой части дифф.уравнения y' = x**3 */
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
                      IntegrationParameters& parameters,
                      unsigned int n)
{
    Table tableB;
    BDF4 bd;
    std::vector<double> x(n, 0.), y(n, 0.);  // инициализация нулями
    RHS rhs;

    double step = parameters.step;

    x[0] = initialState.arg;
    y[0] = initialState.state;   //y0;
    unsigned int srk = sizeof(tableB.table[0])/sizeof(double);   // порядок метода Рунге-Кутты
    unsigned int sbdf = bd.size;                                 // порядок ФДН метода
    std::vector<double>  Ks(srk, 0);

    // разгон (недостающие начальные значения y[i]) методом Рунге-Кутты
    for (unsigned int m = 1; m < sbdf; m++) {
        x[m] = x[0] + m*step;
        y[m] = y[m-1];
        rhs.saa.arg = x[m-1] + tableB.cColumn[0]*step;
        Ks[0] = rhs.calc(rhs.saa);
        for(unsigned int j = 1; j < srk; j++) {
            rhs.saa.arg = x[m-1] + tableB.cColumn[j]*step;
            Ks[j] = rhs.calc(rhs.saa);
        }
        for(unsigned int j = 0; j < sbdf; j++) {
            y[m] = y[m] + Ks[j]*tableB.bString[j]*step;
        }
    }

     // многошаговый метод ФДН
    for (unsigned int m = sbdf; m < n; m++) {
        x[m] = x[0] + m*step;
        rhs.saa.arg = x[m];
        y[m] = bd.alpha[0]*y[m-4] + bd.alpha[1]*y[m-3] + bd.alpha[2]*y[m-2] + bd.alpha[3]*y[m-1] + bd.beta*step*rhs.calc(rhs.saa);

        /* не компилирует undefined reference to BDF4::alpha
        for (unsigned int k = 0; k < sbdf; k++){
            y[m] = y[m] + y[m-sbdf+k]*bd.alpha[k];
        }
        y[m] = y[m] + bd.beta*step*rhs.calc(rhs.saa);
        */
    }

    /* печать решения на экран для проверки
    for (unsigned int i = 0; i < n; i++){
        std::cout << x[i] << ";" << y[i] << ";" << std::endl;
    }
    */

    return y;
}

int main()
{
    double a = 0., b = 5., yInitialPolynom = 0.;     // отрезок на котором находится решение и начальное условие
    std::vector<double> steps {1, 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6};    // шаги интегрирования
    IntegrationParameters params;
    params.epsilon = 1.e-6;
    params.maxIter = 10;

    Polynom pInitial;
    pInitial.saa.arg = a;
    pInitial.saa.state = yInitialPolynom;

    // запись ошибки в файл
    std::ofstream output;
    output.open("bdf4polyErr.csv");
    output << std::setprecision(10);

    // цикл нахождения решения для разных шагов
    for (auto step : steps){
        params.step = step;
        int m = int((b-a)/params.step) + 1;                // количество узлов решения
        std::vector<double> x(m, 0), y(m, 0);           // x,y точного решения
        std::vector<double> y_solution(m, 0);           // вектор вычисленного решения, инициализирован нулями
        for (int i = 0; i < m; i++){
            x[i] = (a + step*i);
            y[i] = func(x[i]);
        }
        y_solution = integrate<RK4Table, Polynom>(pInitial.saa, params, m);         // нахождение вектора решения
        output << params.step << ";" << ErrVectors<double>(y, y_solution) << std::endl;  // записываем ошибку в файл;
    }

    output.close();

    return 0;
}
