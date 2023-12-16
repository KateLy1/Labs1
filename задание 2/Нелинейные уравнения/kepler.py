import pandas as pd
import matplotlib.pyplot as plt
# import math

df = pd.read_csv("EiErr.csv", sep=';', header=None, encoding="cp1251")    # чтение данных из файла
df = df.rename(columns={0: 'ecc', 1: 'i', 2: 'err'})                   # назвали колонки для удобства

fig, ax = plt.subplots()
df.set_index('i', inplace=True)
df.groupby(['ecc'])['err'].plot(legend=True, marker='o')

'''
handles, labels = ax.get_legend_handles_labels()
labels_new = [label.strip('()') for label in labels]
plt.legend(handles, labels_new)
'''

plt.semilogy()
# plt.yscale('log', base=10)
# plt.yticks([1, 1.e-1, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8, 1.e-9, 1.e-10, 1.e-11, 1.e-12, 1.e-13, 1.e-14])

plt.xlabel('Номер итерации i')
plt.ylabel('Логарифм ошибки |Ei - E*|')
plt.title('Зависимость ошибки от номера итерации для разных эксцентриситетов (M=\u03C0/4, точность 10\u207B\u2078)', wrap=True)
plt.show()
plt.close()


# нарисовать степень числа https://ru.stackoverflow.com/questions/1212304/python-Как-нарисовать-степень-числа-используя-функцию-print