#include "integrate.hpp"

int main() {
    const size_t N3 = 3, N5 = 5;  // ïîðÿäîê êâàäðàòóð Ãàóññà äîëæåí áûòü îò 2 äî 5

    std::ofstream output;
    output.open("integralErr.csv");

// std::function<double(double)> sinRef = static_cast<double(*)(double)>(std::sin);

    const double a = 0., b = 10.; // îòðåçîê èíòåãðèðîâàíèÿ [a;b]
    double h = 10, err, integral;

    for (int i = 0; i <= 6; i++) {                                      // öèêë ñ óìåíüøåíèåì øàãà h
        integral = integrate_dx<double, double, N3>(std::sin, a, b, h); // ôóíêöèþ sin(x) ïåðåäàåì â âèäå ïàðàìåòðà
        err = std::abs(std::cos(a) - std::cos(b) - integral);
        output << N3 << ";" << h << ";" << err << std::endl;

        integral = integrate_dx<double, double, N5>(std::sin, a, b, h);
        err = std::abs(std::cos(a) - std::cos(b) - integral);
        output << N5 << ";" << h << ";" << err << std::endl;

        h /= 10;
    }
    output.close();

    return 0;
}

