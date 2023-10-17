import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("newton_interpolation.csv", sep=';', header=None, encoding="cp1251")    # чтение данных из файла
df = df.rename(columns={0: 'N', 1: 'length', 2: 'err', 3: 'Chebyshev'})                   # назвали колонки для удобства

fig, ax = plt.subplots()
df.set_index('length', inplace=True)
df.groupby(['Chebyshev', 'N'])['err'].plot(legend=True, marker='o')

handles, labels = ax.get_legend_handles_labels()
labels_new = [label.strip('()') for label in labels]
plt.legend(handles, labels_new)

plt.xlim((0, 2.1))
plt.xticks([0.0625, 0.125, 0.25, 0.5, 1, 2], ['', '1/8', '1/4', '1/2', '1', '2'])  # надпись 1/16 не входит, поэтому рисуем отдельно
plt.text(-0.06, 0.0000000000011, '1/16')

plt.semilogy()
plt.yticks([1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 1.e-10, 1.e-11])

plt.xlabel('Длина отрезка')
plt.ylabel('Ошибка интерполяции')
plt.title('Зависимость ошибки интерполяции от количества узлов на отрезке. Равномерное и Чебышевское распределение узлов.', wrap=True)
plt.show()
plt.close()
