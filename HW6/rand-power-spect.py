import numpy as np
import matplotlib.pyplot as plt
import time

def linConNumGen(N):
    seed = time.time()
    output = []
    for iii in range(N):
        rand_num = ((22695477 * seed) + 1) % 2**(32)
        seed = rand_num
        output.append(rand_num/2**(32))
    return  np.asarray(output)

def uniformGaussCon(N):
    mu    = 0
    sigma = 1

    u_0 = linConNumGen(N)
    u_1 = linConNumGen(N)
    return np.sqrt(-np.log(u_0)) * np.cos(2*np.pi*u_1), np.sqrt(-np.log(u_0)) * np.sin(2*np.pi*u_1)

def discreteFourier(y):
    size = len(y)
    cs = []
    fs = []
    for iii in range(size):
        temp = 0
        for jjj in range(size):
            temp += y[jjj] * np.exp(-1j * 2 * np.pi * iii * jjj/size)
        cs.append(temp)
        fs.append(2*np.pi*iii/size)
    return np.asarray(cs), np.asarray(fs)

def abs2(x):
    return x.real**2 + x.imag**2

def getRandomWalk(randoms):
    size = len(randoms)
    walk = np.zeros(size)
    for iii in range(1,size):
        walk[iii] = walk[iii-1] + randoms[iii]
    return walk

if __name__ == '__main__':
    outdir = 'outputs/'
    N = 10000
    gauss_1, gauss_2 = uniformGaussCon(N)

    plt.hist(gauss_1, bins=40)
    plt.yscale('log')
    plt.savefig(outdir + 'gauss-hist.png')
    # plt.show()
    plt.clf()

    rand_cs, rand_fs = discreteFourier(gauss_1)
    plt.plot(rand_fs[1:N//2],abs2(rand_cs)[1:N//2])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log k')
    plt.ylabel('log P')
    plt.savefig(outdir + 'gauss-dft.png')
    # plt.show()
    plt.clf()

    walking = getRandomWalk(gauss_1)
    plt.plot(walking)
    plt.savefig(outdir + 'rand-walk.png')
    # plt.show()
    plt.clf()

    walking_cs, walking_fs = discreteFourier(walking)
    plt.plot(walking_fs[1:N//2],abs2(walking_cs)[1:N//2])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log k')
    plt.ylabel('log P')
    plt.savefig(outdir + 'walking-dft.png')
    # plt.show()
