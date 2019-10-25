import numpy as np
import scipy.fftpack as ft
import matplotlib.pyplot as plt
import random as rn
import numpy.fft as ft

def gauss_2d(x,y,sigma,size):
    x = ((x + size//2) % size) - size//2
    y = ((y + size//2) % size) - size//2
    return (1/(sigma*np.sqrt(2*np.pi)))**2 * np.exp(-(.5*(x/sigma)**2) - (.5*(y/sigma)**2))

def readIn(file):
    f=open(file,"r")
    lines=f.readlines()
    result = []
    for x in lines:
        result.append(x.split(' '))
    f.close()
    return np.asarray(result, dtype=float)

img = readIn('blur.txt')
# print(img.shape)
plt.imshow(img, cmap='gray')
plt.savefig('blurr.png')
# plt.show()
plt.clf()
size = img.shape[0]
sigma = 25
xs = range(size)
ys = range(size)
XS, YS = np.meshgrid(xs,ys)
gs = gauss_2d(XS,YS,sigma,size)
plt.imshow(gs, cmap='gray')
plt.savefig('blurrgauss.png')
# plt.show()
plt.clf()

img_ft = ft.rfft2(img)
gs_ft  = ft.rfft2(gs)
print(img_ft.shape, gs_ft.shape)
img_new_ft = np.zeros(img_ft.shape)
for iii in range(size):
    for jjj in range(size//2):
        if ((gs_ft[iii][jjj] < 1e-5) and (gs_ft[iii][jjj] > -1e-5)):
            img_new_ft = img_ft
        else:
            img_new_ft[iii][jjj] = img_ft[iii][jjj]/gs_ft[iii][jjj]
img_new = ft.irfft2(img_new_ft)
print(img_new.shape)
plt.imshow(img_new, cmap='gray')
plt.savefig('blurrenhance.png')
# plt.show()
