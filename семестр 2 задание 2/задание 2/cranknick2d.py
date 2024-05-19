# Реализация схемы Кранка-Николсон - решение уравнения dy/dt = a * d2y/dx2 + f(t,x,y)

import numpy as np
import matplotlib.pyplot as plt
import math


# du/dt = a^2 * (d2u/dx^2 + d2u/dy^2)  u=0 на границах  mxsize - размер матриц A и B  (A*u(n+1) = B*u(n))
def cn_matrix_dirichlet_2d(nn, dx, dt, coeff):
    C0 = coeff * dt / (2 * dx ** 2)
    mxsize = nn * nn
    Amatrix = np.diag(1. + 4. * C0 * np.ones(mxsize)) + \
              np.diag(-C0 * np.ones(mxsize - 1), 1) + \
              np.diag(-C0 * np.ones(mxsize - 1), -1) + \
              np.diag(-C0 * np.ones(mxsize - nn), nn) + \
              np.diag(-C0 * np.ones(mxsize - nn), -nn)

    Bmatrix = np.diag(1. - 4. * C0 * np.ones(mxsize)) + \
              np.diag(C0 * np.ones(mxsize - 1), 1) + \
              np.diag(C0 * np.ones(mxsize - 1), -1) + \
              np.diag(C0 * np.ones(mxsize - nn), nn) + \
              np.diag(C0 * np.ones(mxsize - nn), -nn)

    for i in range(1, mxsize):
        for j in range(1, mxsize):
            if (i % nn == nn - 1 and j % nn == 0) or (i % nn == 0 and j % nn == nn - 1):
                Amatrix[i, j] = 0.
                Bmatrix[i, j] = 0.

    return Amatrix, Bmatrix


# Определение правой части F с учетом ГУ Дирихле
def Fvector(xf, yf, dx, dt, coeff, tf):
    C0 = coeff * dt / (2 * dx ** 2)
    nn = len(xf)
    fv_old = np.zeros((nn, nn))
    fv = np.zeros((nn, nn))
    for i in range(nn):
        for j in range(nn):
            fv_old[i, j] = ((1.+tf*np.pi**2) * xf[j]**2 - 2.*tf) * np.sin(np.pi*yf[i])
            fv[i, j] = ((1.+(tf+dt)*np.pi**2) * xf[j] ** 2 - 2. * (tf+dt)) * np.sin(np.pi*yf[i])
    fv = 0.5 * dt * (fv + fv_old)

    # учет граничных условий
    for i in range(nn):
        fv[i, nn-1] = fv[i, nn-1] + C0 * (tf * np.sin(np.pi * yf[i]) + (tf + dt) * np.sin(np.pi * yf[i]))

    return fv.flatten()    # преобразование двумерного массива в одномерный


# Функция аналитического решения уравнения du/dt = (d2u/dx^2 + d2u/dy^2) +f(t,x,y)
def u_analytical_2d (xf, yf, tf):
    return np.exp(-8.*np.pi**2*tf) * np.sin(2.*np.pi*xf) * np.sin(2.*np.pi*yf) + \
           tf * xf**2 * np.sin(np.pi*yf)


# Сеточная функция аналитического решения уравнения du/dt = (d2u/dx^2 + d2u/dy^2) +f(t,x,y)
def u_analytical_2d_mesh(x, y, t):
    u = np.zeros((len(x), len(y)))
    for i in range(len(x)):           # цикл по всем точкам сетки
        for j in range(len(y)):
            u[i, j] = np.exp(-8.*np.pi**2*t) * np.sin(2.*np.pi*x[j]) * np.sin(2.*np.pi*y[i]) + \
                      t * x[j]**2 * np.sin(np.pi*y[i])
    return u


# Функция ошибки
def max_absolute_error(u, y):
    diff = u - y                            # находим разницу между значениями массивов
    abs_diff = np.absolute(diff)            # находим разницу между значениями массивов по модулю
    max_diff = abs_diff.max()               # находим максимальное значение
    return max_diff


# Итерационный метод Гаусса-Зейделя решения СЛАУ
def gauss_seidel(A, b, x, eps):  # x - начальное приближение
    n = len(A)
    # x = np.zeros_like(b)   # начальное приближение = 0

    while True:
        x_old = np.copy(x)
        for i in range(n):
            temp1 = np.dot(A[i, :i], x[:i])
            temp2 = np.dot(A[i, (i + 1):], x_old[(i + 1):])
            x[i] = (b[i] - temp1 - temp2) / A[i, i]

        if np.sqrt(np.dot(x - x_old, x - x_old)) < eps:
            return x


# Основная программа
def main():
    # Параметры сетки и задачи
    x0, y0 = 0., 0.  # x domain = y domain = [0,1]
    xm, ym = 1., 1.
    Nnodes = [11, 17, 21, 26, 41, 51, 81, 101]    # Количество узлов сетки по x и по y (расчет для разного количества)
    coeff_a = 1.                                  # Коэффициент уравнения
    T = 0.2                                       # Конечное время расчета

    h = np.zeros(len(Nnodes))
    err = np.zeros(len(Nnodes))

    for kn in range(len(Nnodes)):             # Цикл по разным количествам узлов сетки
        h[kn] = (xm - x0) / (Nnodes[kn] - 1)  # Шаг h сетки одинаковый по x и по y с учетом граничных узлов
        nxny = Nnodes[kn] - 2                 # вычитаем граничные узлы
        x = np.linspace(x0+h[kn], xm-h[kn], nxny)       # Узлы равномерной сетки по x без граничных узлов
        y = np.linspace(y0+h[kn], ym-h[kn], nxny)       # Узлы равномерной сетки по y без граничных узлов

        u0 = u_analytical_2d_mesh(x, y, 0.)    # начальные условия u(t=0)
        t = 0                                  # текущее время расчета
        tau = h[kn]                            # Шаг вычисления по времени равен шагу по пространству

        A, B = cn_matrix_dirichlet_2d(nxny, h[kn], tau, coeff_a)  # Заполнение матриц по схеме Кранка-Николсона

        uflat = u0.flatten()    # преобразование двумерного массива решения в одномерный

        # Цикл по времени - получение решения при t = T
        while t < T:
            Fv = Fvector(x, y, h[kn], tau, coeff_a, t)   # правая часть уравнения
            t = round(t + tau, 5)
            # print('t in loop =', t)
            uflat = gauss_seidel(A, np.dot(B, uflat) + Fv, uflat, eps=1.e-3)  # решение СЛАУ методом Гаусса-Зейделя
            # uflat = np.linalg.solve(A, np.dot(B, uflat) + Fv)               # решение СЛАУ не итерационным методом

        u = np.reshape(uflat, (nxny, nxny))               # преобразование одномерного массива решения в двумерный
        uanl = u_analytical_2d_mesh(x, y, t)              # аналитическое решение в узлах сетки
        err[kn] = max_absolute_error(u, uanl)
        print('t=', t, ' h=', h[kn], ' err=', err[kn])

        # Визуализация ошибки решения 3d для каждого h
        fig = plt.figure(figsize=(7, 4))
        ax_3d = fig.add_subplot(projection='3d')
        ax_3d.set_xlabel('x')
        ax_3d.set_ylabel('y')
        ax_3d.set_zlabel('|u_анл - u_числ|')
        ax_3d.set_title("Ошибка при " + " h={:.3f}".format(float(h[kn])))
        xgrid, ygrid = np.meshgrid(x, y)
        # ax_3d.plot_surface(xgrid, ygrid, u, cmap='plasma')
        ax_3d.plot_surface(xgrid, ygrid, np.abs(uanl-u), cmap='plasma')
        plt.show()


    # Визуализация результата
    plt.plot(h, err, '-bo')
    plt.loglog()
    tgangle = round((math.log(err[0]) - math.log(err[2])) / (math.log(h[0]) - math.log(h[2])), 1)
    plt.text(0.07, 0.04, 'k≈' + str(tgangle))
    plt.xlabel('Шаг h')
    plt.ylabel('Ошибка')
    plt.title('Зависимость ошибки от шага')
    plt.show()
    plt.close()


main()

