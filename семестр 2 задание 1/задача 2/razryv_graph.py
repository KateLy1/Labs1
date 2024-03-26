# Отрисовка графиков по данным из файлов
import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image
import sys
import shutil

if len(sys.argv) != 2:
    print("Неверное количество аргументов. \nЗапуск скрипта: python razryv_graph.py KIR или python razryv_graph.py Riman")
    sys.exit(1)

if sys.argv[1] == 'KIR':
    KIR = True
    file_names = ("rho_KIR.txt", "u_KIR.txt", "e_KIR.txt", "p_KIR.txt")   # список файлов с данными для графиков
    dir_name = "./temp_KIR/"                                              # директория для графиков
    scheme = 'Схема КИР'
    yscale_min = (0, 0, 50, 0)
    yscale_max = (14, 250, 180, 11)
elif sys.argv[1] == 'Riman':
    KIR = False
    file_names = ("rho_Riman.txt", "u_Riman.txt", "e_Riman.txt", "p_Riman.txt")  # список файлов с данными для графиков
    dir_name = "./temp_Riman/"                                                   # директория для графиков
    scheme = 'Инвариант Римана'
    yscale_min = (0, 0, 50, 0)
    yscale_max = (14, 600, 180, 11)
else:
    print("Неверный аргумент. Используйте аргумент KIR или Riman для запуска скрипта.")
    sys.exit(1)

# удаление временной директории с файлами, если она существует
if os.path.isdir(dir_name):
    shutil.rmtree(dir_name)
# создание временной директории
os.mkdir(dir_name)

Lx = 10.            # область решения [-L;L] => x=[0;2L] перегородка на расстоянии L
xlabel = 'Метры'
ylabels = ('\u03C1, кг/м3', 'U, м/с', 'E, кДж/кг', 'P, атм')

# чтение файлов, построение графиков во временной директории и создание gif файлов в текущей директории
for j in range(0, len(file_names)):
    with open(file_names[j], "rt") as fcoord:
        nline = 0                               # количество строк в файле
        xlist = fcoord.readline().split(" ")    # пространственная сетка по Х
        x = np.array(xlist, dtype=float)
        x = x - Lx
        frames = []
        while True:
            ylist = fcoord.readline().split(" ")

            if ylist == ['']:
                break

            y = np.array(ylist[1:], dtype=float)   # в y[0] записано время
            frame_n = "Frame_" + str(nline)
            # plt.title(scheme + "   t={:.3f}".format(float(ylist[0])) + ' сек    ' + titles[j])
            plt.title(scheme + "   t={:.3f}".format(float(ylist[0])) + ' сек')
            plt.xlim(-Lx, Lx)
            plt.ylim(yscale_min[j], yscale_max[j])
            plt.xlabel(xlabel)
            plt.ylabel(ylabels[j])
            plt.grid(True)
            plt.plot(x, y)
            file_n = dir_name + frame_n + '_' + file_names[j].rsplit('.', 1)[0] + ".png"  # удалить расширение файла .txt
            plt.savefig(file_n, dpi=200)                 # сохранение изображений в файлы
            plt.close()
            nline += 1

            # Создание gif файлов
            frame = Image.open(file_n)      # Открываем изображение каждого кадра.
            frames.append(frame)            # Добавляем кадр в список с кадрами.
            frames[0].save(                 # Берем первый кадр и в него добавляем оставшееся кадры.
                            file_names[j].rsplit('.', 1)[0] + '.gif',
                            save_all=True,
                            append_images=frames[1:],                   # добавляем со 2-го кадра, т.к.1-й уже есть.
                            optimize=True,                              # сжать палитру, удалив неиспользуемые цвета
                            duration=300,                               # время отображения в миллисекундах
                            loop=0                                      # loop=0 - бесконечный цикл
            )
