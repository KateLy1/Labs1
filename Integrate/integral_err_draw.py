import pandas as pd
import matplotlib.pyplot as plt
import math

df = pd.read_csv("integralErr.csv", sep=';', header=None, encoding="cp1251")    # чтение данных из файла
df = df.rename(columns={0: 'N', 1: 'h', 2: 'err'})                   # назвали колонки для удобства

fig, ax = plt.subplots()
df.set_index('h', inplace=True)
df.groupby(['N'])['err'].plot(legend=True, marker='o')

handles, labels = ax.get_legend_handles_labels()

angle0 = round((math.log(2.13609) - math.log(9.44143e-07))/(math.log(10.) - math.log(1)), 1)
angle1 = round((math.log(0.0223289) - math.log(7.52509e-13))/(math.log(10.) - math.log(1)), 1)

labels_new = [(label.strip('()')) for label in labels]
labels_new[0] = labels_new[0] + ', k≈' + str(angle0)    #'\u00B0'
labels_new[1] = labels_new[1] + ', k≈' + str(angle1)    #'\u00B0'

plt.legend(handles, labels_new)

plt.loglog()

plt.xlabel('Шаг h')
plt.ylabel('Ошибка')
plt.title('Зависимость ошибки вычисления интеграла от шага при N=3 и 5', wrap=True)
plt.show()
plt.close()
