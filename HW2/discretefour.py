import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

def abs2(x):
    return x.real**2 + x.imag**2

def isPower(n):
    n = n / 2
    if n == 2:
        return True
    elif n > 2:
        return isPower(n)
    else:
        return False

def readIn(file):
    f=open(file,"r")
    lines=f.readlines()
    result1=[]
    result2=[]
    for x in lines:
        result1.append(int(x.split('\t')[0]))
        result2.append(float(x.split('\t')[1]))
    f.close()
    return np.asarray(result1), np.asarray(result2)

def discreteFourier(y):
    size = len(y)
    cs = []
    fs = []
    for iii in range(size):
        temp = 0
        for jjj in range(size):
            temp += ys[jjj] * np.exp(-1j * 2 * np.pi * iii * jjj/size)
        cs.append(temp)
        fs.append(2*np.pi*iii/size)
    return np.asarray(cs), np.asarray(fs)

fileName = "sunspots.txt"
xs, ys = readIn(fileName)
size = len(ys)
# wavfile.write('sunspots.wav', 500, ys)

plt.plot(xs, ys)
plt.xlabel('Months')
plt.ylabel('Num. of Sunspots')
plt.savefig('sunspotplot.png')
plt.clf()
# print(24/3100*12)
cs, fs = discreteFourier(ys)

plt.plot(fs[1:size//2],abs2(cs)[1:size//2])
plt.xlabel('freq (1/mo)')
plt.ylabel('DFT of Sunspot |ck|^2')
plt.savefig('sunspotfour.png')
plt.clf()
max_inds = np.argmax(abs2(cs)[1:size//2])
print(fs[max_inds])
# plt.show()
