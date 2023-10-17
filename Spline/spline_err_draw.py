import pandas as pd
import matplotlib.pyplot as plt
import math

df = pd.read_csv("splineErr.csv", sep=';', header=None)    # чтение данных из файла
df = df.rename(columns={0: 'N', 1: 'err'})                   # назвали колонки для удобства

df.plot('N', 'err', legend = None)
plt.loglog()

tgangle = round((math.log(4.18383) - math.log(4275.17))/(math.log(160) - math.log(5)), 1)
plt.text(10, 1600, 'k≈' + str(tgangle))

plt.xlabel('Количество узлов')
plt.ylabel('Ошибка')
plt.title('Сплайн: зависимость ошибки интерполяции от количества узлов', wrap=True)
plt.show()
plt.close()