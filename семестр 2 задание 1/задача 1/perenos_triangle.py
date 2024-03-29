# Реализация схем уголок (параметр 'u' при запуске скрипта)
# и Лакса-Вендрофа (параметр 'lw' при запуске скрипта) для уравнения переноса.
# Периодические граничные условий.
# Запуск скрипта: python perenos.py u или python perenos.py lw

import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Неверное количество аргументов. \nЗапуск скрипта: python perenos.py u или python perenos.py lw")
    sys.exit(1)

if sys.argv[1] == 'lw':
    lw = True
    fig_title =  'Решение уравнения переноса по схеме Лакса-Вендроффа c треугольником '
elif sys.argv[1] == 'u':
    lw = False
    fig_title = 'Решение уравнения переноса по схеме с разностями против потока c треугольником '
else:
    print("Неверный аргумент. Используйте аргумент u или lw для запуска скрипта.")
    sys.exit(1)


# Параметры сетки и задачи
L = 20.0                          # Длина области
Nx = 41                           # Количество точек сетки по пространству
h = L / (Nx - 1)                  # Шаг сетки по пространству  dx=20/40=0.5
x = np.linspace(0, L, Nx)         # Узлы равномерной сетки
v = 1.0                           # Скорость переноса
Kurant = [0.6, 1., 1.01]          # Список чисел Куранта
T = 18.0                          # Время расчета

# Начальные условия y0(x, 0)
# y0 = np.sin(4. * np.pi * x / L)
# Начальные условия - треугольник
y0 = np.copy(x)
y0 = np.where((y0 > 8) & (y0 < 12), y0, 0)
y0 = np.where((y0 > 8) & (y0 <= 10),  0.5 * y0 - 4., y0)
y0 = np.where((y0 > 10) & (y0 < 12), -0.5 * y0 + 6., y0)


# Схема с разностями против потока (уголок) с периодическими граничными условиями с циклом по узлам сетки.
# CFL - число Куранта
def ugolok_periodic(u, nx, CFL):
    unewt = np.zeros_like(u)    # массив для записи решения на следующнм шаге по времени
    for i in range(1, nx):      # вычисление правой границы в узеле nx-1 включено в цикл
        unewt[i] = u[i] - CFL * (u[i] - u[i-1])
    # Периодическое граничное условие - левая граница
    unewt[0] = u[0] - CFL * (u[0] - u[nx - 1])

    return unewt


'''
# Схема Лакса-Вендроффа с периодическими граничными условиями без цикла по узлам сетки
# nx в параметрах для унификации вызова функции
# CFL - число Куранта
def lw_periodic(u, nx, CFL):
    unewt = np.zeros_like(u)  # массив для записи решения на следующнм шаге по времени
    # Центральные точки
    unewt[1:-1] = u[1:-1] - 0.5 * CFL * (u[2:] - u[:-2]) + 0.5 * CFL**2 * (u[2:] - 2*u[1:-1] + u[:-2])
    # Периодическое граничное условие - левая граница
    unewt[0] = u[0] - 0.5 * CFL * (u[1] - u[-1]) + 0.5 * CFL ** 2 * (u[1] - 2 * u[0] + u[-1])
    # Периодическое граничное условие - правая граница
    unewt[-1] = u[-1] - 0.5 * CFL * (u[1] - u[-2]) + 0.5 * CFL ** 2 * (u[1] - 2 * u[-1] + u[-2])
    
    return unewt
'''


# Схема Лакса-Вендроффа с периодическими граничными условиями с циклом по узлам сетки.
# CFL - число Куранта
def lw_periodic(u, nx, CFL):
    unewt = np.zeros_like(u)    # массив для записи решения на следующнм шаге по времени
    for i in range(1, nx - 1):
        unewt[i] = u[i] - 0.5 * CFL * (u[i+1] - u[i-1]) + 0.5 * CFL**2 * (u[i+1] - 2. * u[i] + u[i-1])
    # Периодическое граничное условие - левая граница
    unewt[0] = u[0] - 0.5 * CFL * (u[1] - u[nx - 1]) + 0.5 * CFL**2 * (u[1] - 2. * u[0] + u[nx-1])
    # Периодическое граничное условие - правая граница
    unewt[nx-1] = u[nx-1] - 0.5 * CFL * (u[0] - u[nx-2]) + 0.5 * CFL**2 * (u[0] - 2. * u[nx-1] + u[nx-2])

    return unewt


# Основная программа
figsize = (12, 4)       # размер полотна для рисования графиков
rows = 5                # Таблица графиков 5 х 3
cols = 3
fig, ax = plt.subplots(rows, cols, figsize=figsize, sharex=True, sharey=True)  # таблица графиков с одинаковым масштабом осей
fig.suptitle(fig_title, fontsize=12)
dt_out_results = 5  # промежуток времени между отрисовкой результатов

img_row = 0
img_col = 0

for k in Kurant:            # Цикл по разным числам Куранта
    t = 0                   # текущее время
    y = np.copy(y0)         # начальные условия для каждого нового числа Куранта
    ax[img_row, img_col].set_title('Число Куранта ' + str(k), fontsize=9)
    ax[img_row, img_col].text(16., 0.45, 'T=' + "{:.0f}".format(t), fontsize=8)
    ax[img_row, img_col].plot(x, y)
    time_next_output = 5
    tau = k * h / v         # Шаг вычисления по времени
    nt = int(T / tau) + 1   # Количество временных шагов (с учетом целочисленного преобразования)

    # Цикл по времени
    for j in range(1, nt+1):                # расчет от 1 до nt, чтобы включить конечное время Т
        if lw:
            y = lw_periodic(y, Nx, k)       # решение по схеме Лакса-Вендроффа
        else:
            y = ugolok_periodic(y, Nx, k)   # решение по схеме уголок

        t = j * tau
        # Визуализация результата
        if t >= time_next_output or t >= T:
            time_next_output += dt_out_results
            img_row += 1
            if img_row < rows:
                ax[img_row, img_col].text(16., 0.45, 'T=' + "{:.0f}".format(t), fontsize=8)
                # ax[img_row, img_col].tick_params(axis='both', which='major', labelsize=8)
                ax[img_row, img_col].plot(x, y)

    img_col += 1
    img_row = 0

plt.show()


