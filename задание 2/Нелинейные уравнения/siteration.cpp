#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

template<typename Callable, typename RealType>
const Callable f(RealType x)
{
    return pow(x,2) + pow(std::tan(x),2) - 1.;
}

template<typename Callable, typename RealType>
decltype(auto) solve(
                     const Callable (&func)(RealType),
                     const RealType& tau,
                     const RealType& initialGuess,
                     const unsigned int nIteration
                     )
{
    RealType x0, x = initialGuess, tol = 1.e-6;
    unsigned int iter = 0;

    do {
        x0 = x;
        x = x0 + tau*f<Callable>(x0);   // МПИ с релаксацией
        iter++;

        if (iter > nIteration) {
            throw nIteration;
        }
        std::cout << "Iteration " << iter << ": " << "x=" << x << "   ";
        std::cout << "|x(i+1) - x(i)| = " << std::abs(x - x0) << std::endl;
    } while (std::abs(x - x0) > tol);

    return x;
}

int main()
{
    double x, tau = -0.15, initialGuess = 0.5;      // 1-й корень (положительный)
    //double x, tau = 0.15, initialGuess = -0.5;    // 2-й корень (отрицательный)
    unsigned int nIteration = 20  ;
    std::cout << std::endl << std::setprecision(7);

    try {
        x = solve<double>(f, tau, initialGuess, nIteration);
        std::cout << "Solution x = " << x;

        }
        catch (unsigned int& nIter) {           // исключение, если решение не нашлось за nIteration
            std::cout << "For initialGuess=" << initialGuess << " solution was not found for " << nIter << " iterations" << std::endl;
        }

    return 0;

}


