import numpy as np
import random as rn
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def createMagArray(M):
    magarray = np.zeros((M,M))
    for iii in range(M):
        for jjj in range(M):
            if rn.uniform(0,1) < 0.5:
                magarray[iii,jjj] = 1
            else:
                magarray[iii,jjj] = -1
    return magarray

def calcEnergy(magarray,J):
    size = len(magarray)
    energy = 0
    for ggg in range(size//2):
        for hhh in range(size):
            if hhh+1 < size:
                energy += (- J * magarray[2*ggg,hhh] * magarray[2*ggg,hhh+1])
            if ggg+1 < size:
                energy += (- J * magarray[2*ggg,hhh] * magarray[(2*ggg)+1,hhh])
            if ggg-1 > -1:
                energy += (- J * magarray[2*ggg,hhh] * magarray[(2*ggg)-1,hhh])
    return energy

def differenceEnergy(diffcell_i,diffcell_j,magarray,J):
    size = len(magarray)
    iii = diffcell_i
    jjj = diffcell_j
    energy = 0
    if (iii+1 < size):
        energy += (-J * -1 * magarray[iii,jjj] * magarray[iii+1,jjj])
    if (iii-1 > -1):
        energy += (-J * -1 * magarray[iii,jjj] * magarray[iii-1,jjj])
    if (jjj+1 < size):
        energy += (-J * -1 * magarray[iii,jjj] * magarray[iii,jjj+1])
    if (jjj-1 > 0):
        energy += (-J * -1 * magarray[iii,jjj] * magarray[iii,jjj-1])
    return energy

def getTotalMag(array):
    return np.sum(array)

plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

FFMpegWriter = ani.writers['ffmpeg']
metadata = dict(title='isingmodel', artist='MattDaunt',
                comment='ising model of magnetic array of particles')
writer = FFMpegWriter(fps=40, metadata=metadata)
outdir = 'outputs/'
fig = plt.figure()

M = 100
J = 1
Temp = 3
magarray = createMagArray(M)
Etotal = calcEnergy(magarray,J)
max_iter = 1000000

total_mag = np.zeros(max_iter)

temp_array = [Temp * np.exp(-4 * iii / max_iter) for iii in range(max_iter)]
# THE ISING MODEL
with writer.saving(fig, outdir + 'ising_model_t' + str(Temp) + '.mp4', 100):
    for iter in range(max_iter):
        total_mag[iter] = getTotalMag(magarray)

        new_cell_i = rn.randint(0,M-1)
        new_cell_j = rn.randint(0,M-1)
        Ediff = differenceEnergy(new_cell_i,new_cell_j,magarray,J)
        # print(Ediff)
        if Ediff < 0:
            # print('smol E')
            Etotal += Ediff
            magarray[new_cell_i,new_cell_j] *= -1
            if (iter % 1000) == 0:
                plt.imshow(magarray)
                writer.grab_frame()
                plt.clf()
        elif rn.uniform(0,1) < np.exp(-Ediff/Temp):
            # print('probably added')
            Etotal += Ediff
            magarray[new_cell_i,new_cell_j] *= -1
            if (iter % 1000) == 0:
                plt.imshow(magarray)
                writer.grab_frame()
                plt.clf()
        # else:
            # print('rejected')
plt.plot(total_mag)
plt.ylabel('Net Magnetization')
plt.xlabel('iterations')
plt.savefig(outdir + 'ising-netmag-t' + str(Temp) + '.png')
plt.show()
# plt.imshow(magarray)
# plt.show()
