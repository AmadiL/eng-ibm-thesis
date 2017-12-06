import numpy as np
import os

with open('testy/test2/raw.csv') as f:
    sample = np.fromstring(f.read(), dtype=np.float64, sep=os.linesep)


fs = sample[0]
y = sample[1:]
N = len(y)
t = np.arange(N)/fs
noise = np.random.normal(0, 0.01, N)
s = {}
for i in range(1, 11):
    s[i] = y + np.random.normal(0, 0.1*i, N)
    with open('testy/test2/gaussian{}.csv'.format(i), 'w') as f:
        f.write(str(fs)+'\n')
        f.writelines(list(map(lambda x: str(x)+'\n', s[i])))
