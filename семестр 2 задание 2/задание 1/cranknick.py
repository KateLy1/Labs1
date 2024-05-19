# Реализация схемы Кранка-Николсон решение уравнения dy/dt = a * d2y/dx2 + f(t,x)

import numpy as np
import matplotlib.pyplot as plt
import math


# du/dt = a^2 * d2u/dx^2  u(0,t) = u(L,t) = 0  u(0) = u0  mxsize - размер матриц A и B  (A*u(n+1) = B*u(n))
def cn_matrix_dirichlet(mxsize, dx, dt, coeff):
    C0 = coeff * dt / (2 * dx ** 2)
    Amatrix = np.diag(1 + 2 * C0 * np.ones(mxsize)) + \
              np.diag(-C0 * np.ones(mxsize-1), 1) + \
              np.diag(-C0 * np.ones(mxsize-1), -1)
    Bmatrix = np.diag(1 - 2 * C0 * np.ones(mxsize)) + \
              np.diag(C0 * np.ones(mxsize-1), 1) + \
              np.diag(C0 * np.ones(mxsize-1), -1)

    return Amatrix, Bmatrix


# Функция аналитического решения уравнения du/dt = a^2 * d2u/dx^2  u(x0,t) = u(xM,t) = 0  u(t=0) = u0 = const
# Решение существует в виде суммы бесконечного ряда.
# eps - точность аналитического решения
def u_analytical(x, t, u0, a, eps):
    u = np.zeros_like(u0)
    for i in range(len(x)):           # цикл по всем точкам сетки
        k = 0
        term_k = 2. * eps             # фиктивное значения для запуска цикла while
        while np.abs(term_k) > eps:   # term_k - k-й член ряда
            term_k = (4. * u0[i] / np.pi) * (1. / (2 * k + 1)) * np.exp(-(a * (2 * k + 1) * np.pi / x[-1]) ** 2 * t) * \
                     np.sin((2 * k + 1) * np.pi * x[i] / x[-1])
            u[i] += term_k
            k += 1
    return u


# Функция ошибки
def max_absolute_error(u, y):
    diff = u - y                            # находим разницу между значениями массивов
    abs_diff = np.absolute(diff)            # находим разницу между значениями массивов по модулю
    max_diff = abs_diff.max()               # находим максимальное значение
    return max_diff


# Основная программа
def main():
    # Параметры сетки и задачи
    x0 = 0.                                        # Начало области [x0; xM]
    xM = 10.                                       # Конец области [x0; xM]
    Nnodes = [6, 11, 21, 41, 101, 201, 401, 1001]  # Количество узлов сетки по x (расчет для разного количества узлов)
    coeff_a = 1.                                  # Коэффициент уравнения (к-нт диффузии для уравнения теплопроводности)
    T = 6.                                        # Конечное время расчета

    h = np.zeros(len(Nnodes))
    err = np.zeros(len(Nnodes))

    for kn in range(len(Nnodes)):            # Цикл по разным количествам узлов сетки
        Nx = Nnodes[kn]
        h[kn] = (xM - x0) / (Nx - 1)         # Шаг h сетки по пространству
        x = np.linspace(x0, xM, Nx)          # Узлы равномерной сетки
        u0 = np.full(Nx, 5.)                 # начальные условия u(t=0) = 5
        boundary_dirichlet = np.full(2, 0.)  # Граничные условия Дирихле u(x0,t) = u(xM,t) = 0.
        u0[0] = boundary_dirichlet[0]
        u0[-1] = boundary_dirichlet[1]
        u = np.copy(u0)                      # начальные условия (t=0)
        t = 0                                # текущее время
        tau = h[kn]                          # Шаг вычисления по времени

        A, B = cn_matrix_dirichlet(len(u)-2, h[kn], tau, coeff_a)  # Заполнение матриц по схеме Кранка-Николсона

        # Цикл по времени - получение решения при t = T
        while t < T:
            t = t + tau
            u[1:-1] = np.linalg.solve(A, np.dot(B, u[1:-1]))  # решение СЛАУ в центральных точках

        y = u_analytical(x, t, u0, coeff_a, eps=1.e-12)                   # аналитическое решение в узлах сетки
        err[kn] = max_absolute_error(u, y)
        print('h=', h[kn], 'err=', err[kn])

    # Визуализация результата
    plt.plot(h, err, '-bo')
    plt.loglog()
    tgangle = round((math.log(err[1]) - math.log(err[3])) / (math.log(h[1]) - math.log(h[3])), 1)
    plt.text(0.4, 0.0015, 'k≈' + str(tgangle))
    plt.xlabel('Шаг h')
    plt.ylabel('Ошибка')
    plt.title('Зависимость ошибки от шага')
    plt.show()
    plt.close()


main()

