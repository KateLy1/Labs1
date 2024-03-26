# Задача Римана о распаде разрыва
# Решение с помощью схемы Куранта-Изаксона-Риса (КИР) или
# с использованием подхода с переходом к инвариантам Римана
# в зависимости от аргументов скрипта
# В текущей директории создаются текстовые файлы с записями решения каждую миллисекунду в интервале времени от 0 до Т
# Запуск скрипта: python razryv.py KIR или python razryv.py Riman

import numpy as np
import math
import sys

if len(sys.argv) != 2:
    print("Неверное количество аргументов. \nЗапуск скрипта: python razryv.py KIR или python razryv.py Riman")
    sys.exit(1)

if sys.argv[1] == 'KIR':
    KIR = True
    file_names = ("rho_KIR.txt", "u_KIR.txt", "e_KIR.txt", "p_KIR.txt")          # список файлов для записи решения КИР
elif sys.argv[1] == 'Riman':
    file_names = ("rho_Riman.txt", "u_Riman.txt", "e_Riman.txt", "p_Riman.txt")  # список файлов для записи решения Риман
    KIR = False
else:
    print("Неверный аргумент. Используйте аргумент KIR или Riman для запуска скрипта.")
    sys.exit(1)

Lx = 10.            # область решения [-L;L] => x=[0;2L] перегородка на расстоянии L
T = 0.02            # решение от 0 до Т сек
t = 0.              # начальное время = 0
Nx_half = 50
Nx = 2 * Nx_half   # количество узлов должно быть четным (значение параметров не определено на перегородке, поэтому узел не должен попадать на перегородку)
gamma = 5./3.

# пространственная сетка по Х с узлами l
xl = np.linspace(0, 2.*Lx, Nx)       # равномерная сетка по Х
h = xl[1] - xl[0]		            # шаг по координате Х


def matrixA(u, e):
    return np.array([[0., 1, 0.], [-u**2, 2*u, gamma-1], [-gamma*u*e, gamma*e, u]])


def OmegaT(u, e):
    csq = gamma * (gamma - 1) * e
    c = math.sqrt(csq)
    return np.array([[-u*c, c, gamma-1], [-csq, 0., gamma-1], [u*c, -c, gamma-1]])


def inverseOmegaT(u, e):
    arr = OmegaT(u, e)
    return np.linalg.inv(arr)


def absLamda(u, e):
    csq = gamma * (gamma - 1) * e
    c = math.sqrt(csq)
    return np.array([[math.fabs(u+c), 0., 0.], [0., math.fabs(u), 0.], [0., 0., math.fabs(u-c)]])


def Lamda(u, e):
    csq = gamma * (gamma - 1) * e
    c = math.sqrt(csq)
    return np.array([[u+c, 0., 0.], [0., u, 0.], [0., 0., u-c]])


# начальные условия состояния газа слева и справа от перегородки
# для консервативных переменных w(rho, rho*u, rho*e)
def initial_conditions():
    # слева от перегородки
    rho_left = 13.
    p_left = 101325.011 * 10.  # 10 атмосфер выраженных в Паскалях
    u_left = 0.
    
    # справа от перегородки
    rho_right = 1.3
    p_right = 101325.011 * 1.  # 1 атмосфера выраженная в Паскалях
    u_right = 0.
    
    wl = np.zeros((3, Nx), dtype="float64")
    wl[0][:Nx_half] = rho_left
    wl[1][:Nx_half] = rho_left * u_left
    wl[2][:Nx_half] = p_left / (gamma - 1)

    wl[0][Nx_half:] = rho_right
    wl[1][Nx_half:] = rho_right * u_right
    wl[2][Nx_half:] = p_right / (gamma - 1)     # p/(gamma-1) = rho*e

    return wl


# Схема КИР для вектора консервативных переменных w(rho, rho*u, rho*e)
# dt - шаг по времени (в общем случае неравномерный и зависит от числа Куранта для выполнения условия устойчивости схемы)
def calc_KIR(wl, dt):
    max_lmd = 0.
    wl1 = np.copy(wl)           # создание копии текущего вектора для вектора на следующем временном шаге
    for i in range(1, Nx-1):    # цикл по пространственной сетке
        ul = wl[1][i] / wl[0][i]
        el = wl[2][i] / wl[0][i]
        Al = matrixA(ul, el)

        wldi[0] = (wl[0][i+1] - wl[0][i-1]) / (2.*h)
        wldi[1] = (wl[1][i+1] - wl[1][i-1]) / (2.*h)
        wldi[2] = (wl[2][i+1] - wl[2][i-1]) / (2.*h)
        Ali = Al @ wldi

        OmTl = OmegaT(ul, el)
        Laml = absLamda(ul, el)
        max_lmd = max(max_lmd, max([max(m) for m in Laml]))     # нахождение макс. числа матрицы
        invOmTl = inverseOmegaT(ul, el)
        wlddi[0] = (wln[0][i+1] - 2.*wln[0][i] + wln[0][i-1]) / (2.*h)
        wlddi[1] = (wln[1][i+1] - 2.*wln[1][i] + wln[1][i-1]) / (2.*h)
        wlddi[2] = (wln[2][i+1] - 2.*wln[2][i] + wln[2][i-1]) / (2.*h)
        Oli = (invOmTl @ Laml @ OmTl) @ wlddi

        # основное вычисление вектора w(rho, rho*u, rho*e)
        wl1[:, i] = wl[:, i] - dt * Ali[:] + dt * Oli[:]

    # граничные условия grad = 0
    wl1[0][0] = wl1[0][1]
    wl1[1][0] = wl1[1][1]
    wl1[2][0] = wl1[2][1]
    wl1[0][Nx-1] = wl1[0][Nx-2]
    wl1[1][Nx-1] = wl1[1][Nx-2]
    wl1[2][Nx-1] = wl1[2][Nx-2]

    return wl1, max_lmd


# Использование подхода с переходом к инвариантам Римана для вектора консервативных переменных w(rho, rho*u, rho*e)
# dt - шаг по времени (в общем случае неравномерный и зависит от числа Куранта для выполнения условия устойчивости схемы)
def calc_invariant_Riman(wl, dt):
    max_lmd = 0.
    wl1 = np.copy(wl)           # создание копии текущего вектора для вектора на следующем временном шаге
    for i in range(1, Nx-1):    # цикл по пространственной сетке
        ul = wl[1][i] / wl[0][i]
        el = wl[2][i] / wl[0][i]
        ulm1 = wl[1][i-1] / wl[0][i-1]  # u[l-1]
        elm1 = wl[2][i-1] / wl[0][i-1]
        ulp1 = wl[1][i+1] / wl[0][i+1]  # u[l+1]
        elp1 = wl[2][i+1] / wl[0][i+1]
        OmTl = OmegaT(ul, el)
        invOmTl = inverseOmegaT(ul, el)
        OmTlm1 = OmegaT(ulm1, elm1)
        OmTlp1 = OmegaT(ulp1, elp1)
        Laml_abs = absLamda(ul, el)
        Laml = Lamda(ul, el)

        wldi = (OmTl @ wl[:, i] - OmTlm1 @ wl[:, i-1]) / (2. * h)
        Lm1 = invOmTl @ (Laml + Laml_abs) @ wldi

        wlddi = (OmTlp1 @ wl[:, i+1] - OmTl @ wl[:, i]) / (2. * h)
        Lm2 = invOmTl @ (Laml - Laml_abs) @ wlddi

        max_lmd = max(max_lmd, max([max(m) for m in Laml_abs]))     # нахождение макс. числа матрицы

        # основное вычисление вектора w(rho, rho*u, rho*e)
        wl1[:, i] = wl[:, i] - dt * Lm1[:] - dt * Lm2[:]

    # граничные условия grad = 0
    wl1[0][0] = wl1[0][1]
    wl1[1][0] = wl1[1][1]
    wl1[2][0] = wl1[2][1]
    wl1[0][Nx-1] = wl1[0][Nx-2]
    wl1[1][Nx-1] = wl1[1][Nx-2]
    wl1[2][Nx-1] = wl1[2][Nx-2]

    return wl1, max_lmd


# Компоненты вектора wln: wln[0]=rho, wln[1]=rho*u, wln[2]=rho*e, p=(gamma-1)*wln[2]
wln0 = initial_conditions()             # начальные условия вектора w(rho, rho*u, rho*e)
p = wln0[2] * (gamma-1) / 101325.011   # давление p=(gamma-1)*rho*e в атмосферах
wln = np.copy(wln0)                     # w на текущем временном шаге
wln1 = np.copy(wln)                     # w на следующем временном шаге
wldi = np.zeros(3, dtype="float64")     # w для промежуточных вычислений
wlddi = np.zeros(3, dtype="float64")    # w для промежуточных вычислений

# открытие файлов в текущей директории для записи решения по каждой компоненте вектора w(3) и давления p
fid = []    # список id файлов
for fn in file_names:
    fid.append(open(fn, "w", newline=""))  # файл создается или перезаписывается, если существует

# запись пространственной сетки и начальных условий в файлы для каждой компоненты вектора w и давления
for i in range(0, len(fid)):
    xl.tofile(fid[i], sep=' ')
    fid[i].write('\n')
    fid[i].write(str(t)+" ")          # значения t=0 (начальные условия)
    if i < 2:                               # rho и u
        wln0[i].tofile(fid[i], sep=' ')     # начальные условия
        fid[i].write('\n')
    elif i == 2:                            # e - энергия
        wln0[2] = wln0[2] / wln0[0]         # rho*e/rho
        wln0[2] = wln0[2] / 1000.           # кДж
        wln0[i].tofile(fid[i], sep=' ')     # начальные условия
        fid[i].write('\n')
    elif i == 3:                        # для давления
        p.tofile(fid[i], sep=' ')      # начальные условия
        fid[i].write('\n')


# Вычисление решения на каждом шаге по времени
Kurant_max = 0.01  # Макс.число Куранта для этой задачи. Kurant = tau*max|lamda(i)|/h < 1 (для нелинейных <<1)
tau_init = 1.e-6
tau = tau_init
pwr = 0

# цикл по времени
while t <= T:
    t = round(t + tau, 7)

    # вычисление вектора w на следующем временном шаге и макс.собственного числа для схемы КИР или инварианта Римана
    if KIR:
        wln1, max_lamda = calc_KIR(wln, tau)
    else:
        wln1, max_lamda = calc_invariant_Riman(wln, tau)

    # Контроль условия устойчивости.
    Kurant_num = tau * max_lamda / h
    if Kurant_num > Kurant_max:             # Для данной нелинейной задачи Kurant должен не превышать Kurant_max
        pwr = int(math.log10(Kurant_num / Kurant_max)) + 1    # на сколько порядков число Куранта выше макс. Куранта
        tau = tau / (10**pwr)                              # уменьшаем шаг по времени на нужное число порядков

    wln = np.copy(wln1)
    p = wln[2] * (gamma-1) / 101325.011   # давление в атмосферах

    # запись решения на шаге по времени в файлы для каждой компоненты вектора w
    if int(t*1.e+6) % 1000 == 0:           # запись каждую миллисекунду
        print('t=', t*1000, 'мсек', '  Kurant=', "{:6f}".format(Kurant_num), '  tau=', tau, '  pwr=', pwr)

        for i in range(0, len(fid)):
            fid[i].write(str("{:.6f}".format(t)) + " ")  # значение времени
            if i == 0:  # запись rho - плотность
                wln[i].tofile(fid[i], sep=' ')
                fid[i].write('\n')
            elif i == 1:  # запись u - скорость
                wln1[i] = wln[i] / wln[0]          # временное использование переменной wln1 для вычисления консервативной переменной u
                wln1[i].tofile(fid[i], sep=' ')
                fid[i].write('\n')
            elif i == 2:  # запись e - энергия
                wln1[i] = wln[i] / wln[0]        # временное использование переменной wln1 для вычисления консервативной переменной e
                wln1[i] = wln1[i] / 1000.        # в кДж
                wln1[i].tofile(fid[i], sep=' ')
                fid[i].write('\n')
            elif i == 3:  # запись p - давление
                p.tofile(fid[i], sep=' ')
                fid[i].write('\n')

# закрытие файлов
for i in range(0, len(fid)):
    fid[i].close()



