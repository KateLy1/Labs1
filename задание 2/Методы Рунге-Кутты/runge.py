import pandas as pd
import matplotlib.pyplot as plt
import math

df1 = pd.read_csv("polynomErr.csv", sep=';', header=None)    # чтение данных из файла
df1 = df1.rename(columns={0: 'step', 1: 'pol'})                   # назвали колонки для удобства
#df.set_index('step', inplace=True)

df2 = pd.read_csv("oscErr.csv", sep=';', header=None)
df2 = df2.rename(columns={0: 'step', 1: 'osc'})
#df2.set_index('step', inplace=True)

#ax = df1.plot(style=['b','y','g'])
# df2.plot(ax=ax, style=['b','y','g'], linestyle='--')

ax = df1.plot('step', 'pol', legend = None)
df2.plot('step', 'osc', legend = None, ax=ax)
plt.loglog()
plt.text(2.6e-2, 1.e-5, 'y\'\' = -y')
plt.text(2.5e-6, 2.2e-11, 'y\' = x\u00B3')

tgangle1 = round((math.log(7.105427358e-15) - math.log(1.261923899e-11))/(math.log(1) - math.log(1e-6)), 1)
plt.text(5e-5, 1e-11, 'k≈' + str(tgangle1))
tgangle2 = round((math.log(0.03269948605) - math.log(2.681188604e-14))/(math.log(1) - math.log(1e-03)), 1)
plt.text(1e-2, 1.18e-7, 'k≈' + str(tgangle2))


plt.xlabel('Шаг')
plt.ylabel('Ошибка')
plt.title('Метод Рунге-Кутты 4-го порядка: зависимость ошибки от шага', wrap=True)
plt.show()
plt.close()


# нарисовать степень числа https://ru.stackoverflow.com/questions/1212304/python-Как-нарисовать-степень-числа-используя-функцию-print