//https://en.wikipedia.org/wiki/Kepler%27s_equation
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol)
{
    double E0, E = meanAnomaly;     // начальное приближение записано в переменную E
    double f, df;                   // функция f(E)=0 и ее производная
    unsigned int iter = 0;
    std::vector<double> Ei(maxIter+1, 100.);  // Ei - вектор для записи значения Е на каждой итерации (заполнен значением 100)

    Ei[0] = E;                                // первый элемент - начальное приближение (нулевая итерация)

    do {
        E0 = E;
        f = E0 - ecc*std::sin(E0) - meanAnomaly;
        df = 1. - ecc*std::cos(E0);
        E = E0 - f/df;               // метод Ньютона (касательных)
        iter++;

        if (iter > maxIter) {
            throw maxIter;
        }

        Ei[iter] = E;

    } while (std::abs(E - E0) > tol);

    // запись в файл разности полученного решения Е и значений на каждой итерации Ei (без последнего значения ошибки, которое равно 0)
    for (unsigned j = 0; j < iter; j++) {
        std::cout << ecc << ";" << j <<";" << std::abs(E - Ei[j]) << ";" << std::endl;
    }
    // std::cout << E << std::endl;
    return E;
}

int main()
{
    std::vector<double> ecc {0.1, 0.2, 0.5, 0.8};   // эксцентриситеты
    double meanAnomaly = M_PI/4., tol = 1.e-8;
    unsigned int maxIter = 10;

    // организация потока вывода cout в файл csv
    std::ofstream out("EiErr.csv");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf());                //redirect std::cout to file
    std::ofstream output;

    for (unsigned int k = 0; k < ecc.size(); k++){
        try {
            keplerSolver(ecc[k], meanAnomaly, maxIter, tol);
        }
        catch (unsigned int& nIter) {           // исключение, если решение не нашлось - запись об этом в файл и продолжение программы
            std::cout << "For ecc=" << ecc[k] << " solution was not found for " << nIter << " iterations" << std::endl;
        }
    }

    std::cout.rdbuf(coutbuf); //reset to standard output again

    return 0;


}
