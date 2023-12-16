import pandas as pd
import matplotlib.pyplot as plt
import math

df1 = pd.read_csv("bdf4polyErr.csv", sep=';', header=None)    # чтение данных из файла
df1 = df1.rename(columns={0: 'step', 1: 'pol'})                   # назвали колонки для удобства
#df.set_index('step', inplace=True)

df2 = pd.read_csv("bdf4oscErr.csv", sep=';', header=None)
df2 = df2.rename(columns={0: 'step', 1: 'osc'})

#ax = df1.plot(style=['b','y','g'])
# df2.plot(ax=ax, style=['b','y','g'], linestyle='--')

ax = df1.plot('step', 'pol', legend = None)
df2.plot('step', 'osc', legend = None, ax=ax)
plt.loglog()
plt.text(1.1e-2, 1.e-5, 'y\'\' = -y')
plt.text(2.5e-1, 2.2e-13, 'y\' = x\u00B3')

tgangle1 = round((math.log(7.105427358e-15) - math.log(1.672049166e-10))/(math.log(1) - math.log(1e-06)), 1)
plt.text(1.18e-2, 1.18e-12, 'k≈' + str(tgangle1))
tgangle2 = round((math.log(0.1603385243) - math.log(6.967759703e-13))/(math.log(1) - math.log(1e-03)), 1)
plt.text(5e-3, 1.18e-7, 'k≈' + str(tgangle2))

plt.xlabel('Шаг')
plt.ylabel('Ошибка')
plt.title('Неявный ФДН 4-го порядка: зависимость ошибки от шага', wrap=True)
plt.show()
plt.close()

# нарисовать степень числа https://ru.stackoverflow.com/questions/1212304/python-Как-нарисовать-степень-числа-используя-функцию-print