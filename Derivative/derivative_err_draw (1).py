import pandas as pd
import matplotlib.pyplot as plt
import math

df = pd.read_csv("derivativeErr.csv", sep=';', header=None, encoding="cp1251")    # чтение данных из файла
df = df.rename(columns={0: 'N', 1: 'h', 2: 'err'})                   # назвали колонки для удобства

fig, ax = plt.subplots()
df.set_index('h', inplace=True)
df.groupby(['N'])['err'].plot(legend=True, marker='o')

handles, labels = ax.get_legend_handles_labels()
labels_new = [label.strip('()') for label in labels]

tgangle0 = round((math.log(0.369289) - math.log(2.27433e-07))/(math.log(1.) - math.log(0.01)), 1)
tgangle1 = round((math.log(0.357145) - math.log(1.37052e-09))/(math.log(1.) - math.log(0.01)), 1)
tgangle2 = round((math.log(0.0816284) - math.log(4.5528e-12))/(math.log(1.) - math.log(0.01)), 1)

labels_new = [(label.strip('()')) for label in labels]
labels_new[0] = labels_new[0] + ', k≈' + str(tgangle0)    #'\u00B0'
labels_new[1] = labels_new[1] + ', k≈' + str(tgangle1)    #'\u00B0'
labels_new[2] = labels_new[2] + ', k≈' + str(tgangle2)    #'\u00B0'

plt.legend(handles, labels_new)

plt.loglog()

plt.xlabel('Шаг h')
plt.ylabel('Ошибка')
plt.title('Зависимость ошибки вычисления производной от шага при N=3,4,5', wrap=True)
plt.show()
plt.close()
