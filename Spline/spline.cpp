#include "spline.hpp"

int main() {
	double a = 0, b = 10;                               // отрезок [a;b]
    std::vector<int> N_nodes {5, 10, 20, 40, 80, 160};  // различное количество узлов сплайна
    std::vector<double> x, y;                           // x, y -координаты узлов
	int np = 1000;                                      // количество точек интерпол€ции
	std::vector<double> xp, yp, y_func;                 // координаты точек интерпол€ции и значени€ непрерывной функции в xp
	std::vector<double> splineErr;                      // вектор дл€ записи ошибки интерпол€ции при разном количестве узлов
    double s = (b-a)/(np-1);
    for (int i = 0; i < np; i++) {
        xp.push_back(a + s*i);
        y_func.push_back(std::exp(a + s*i));            // значени€ функции в точках интерпол€ции
    }

    // цикл построени€ сплайна по разному количеству узлов
    for (auto n : N_nodes){
        double h = (b-a)/(n-1);   // разбиение отрезка
        for (int i = 0; i < n; i++){
            x.push_back(a + h*i);
            y.push_back(std::exp(a + h*i));
        }
        CubicSpline<double, double> spline(x, y);
        yp = spline.interpolate(xp);
        splineErr.push_back(ErrVectors<double>(yp, y_func));  // записываем ошибку интерпол€ции в вектор;

        std::vector<double>().swap(x);                        // удал€ем все элементы вектора x и освобождаем пам€ть дл€ инициализации с другим количеством узлов
        std::vector<double>().swap(y);                        // удал€ем все элементы вектора y и освобождаем пам€ть дл€ инициализации с другим количеством узлов
    }

    // запись результата в файл
    std::ofstream output;
    output.open("splineErr.csv");
    for (unsigned int i = 0; i < N_nodes.size(); i++) {
        output << N_nodes[i] << ";" << splineErr[i] << std::endl;
    }
    output.close();

    return 0;

}

