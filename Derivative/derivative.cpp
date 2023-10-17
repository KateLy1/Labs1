#include "derivative.hpp"

int main() {

    std::vector<double> points3 {-1, 1, 2}, points4 {-1, 1, 2, 3}, points5 {-2, -1, 1, 2, 3};
    DerivativeCoef<double> dcoef3, dcoef4, dcoef5;

    dcoef3 = calcDerivativeCoef<double>(points3);
    printCoefs(dcoef3, points3);
    dcoef4 = calcDerivativeCoef<double>(points4);
    printCoefs(dcoef4, points4);
    dcoef5 = calcDerivativeCoef<double>(points5);
    printCoefs(dcoef5, points5);


// запись результата в файл
    std::ofstream output;
    output.open("derivativeErr.csv");

    double dt, err, x0 = 1, h = 1;
    for (int i = 0; i <= 15; i++) {
        dt = calcDerivative1st<double>(std::exp, dcoef3, points3, x0, h); // exp передаем в виде параметра
        err = std::abs(std::exp(x0) - dt);
        output << points3.size() << ";" << h << ";" << err << std::endl;

        dt = calcDerivative1st<double>(std::exp, dcoef4, points4, x0, h);
        err = std::abs(std::exp(x0) - dt);
        output << points4.size() << ";" << h << ";" << err << std::endl;

        dt = calcDerivative1st<double>(std::exp, dcoef5, points5, x0, h);
        err = std::abs(std::exp(x0) - dt);
        output << points5.size() << ";" << h << ";" << err << std::endl;

        h /= 10;
    }

    return 0;
}
