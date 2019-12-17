import numpy as np
import random as rn
import matplotlib.pyplot as plt
import matplotlib.animation as ani

def calcEnergy(latice):
    size = len(latice)
    energy = 0
    for iii in range(size):
        for jjj in range(size):
            if latice[iii,jjj] != 0:
                energy -= 1
    energy /= 2
    return energy

def getNewLatice(latice,site_i,site_j,adj_i,adj_j,iter):
    new_latice = np.array(latice)
    size       = len(latice)
    if new_latice[site_i,site_j] != 0:
        dimer_1 = new_latice[site_i,site_j]
        new_latice[site_i,site_j] = 0
        if site_i+1 < size and new_latice[site_i+1,site_j] == dimer_1:
            new_latice[site_i+1,site_j] = 0
        elif site_i > 0 and new_latice[site_i-1,site_j] == dimer_1:
            new_latice[site_i-1,site_j] = 0
        elif site_j+1 < size and new_latice[site_i,site_j+1] == dimer_1:
            new_latice[site_i,site_j+1] = 0
        elif site_j > 0 and new_latice[site_i,site_j-1] == dimer_1:
            new_latice[site_i,site_j-1] = 0
    if new_latice[adj_i,adj_j] != 0:
        dimer_2 = new_latice[adj_i,adj_j]
        new_latice[adj_i,adj_j] = 0
        if adj_i+1 < size and new_latice[adj_i+1,adj_j] == dimer_2:
            new_latice[adj_i+1,adj_j] = 0
        elif adj_i > 0 and new_latice[adj_i-1,adj_j] == dimer_2:
            new_latice[adj_i-1,adj_j] = 0
        elif adj_j+1 < size and new_latice[adj_i,adj_j+1] == dimer_2:
            new_latice[adj_i,adj_j+1] = 0
        elif adj_j > 0 and new_latice[adj_i,adj_j-1] == dimer_2:
            new_latice[adj_i,adj_j-1] = 0
    iter += 1
    new_latice[site_i,site_j] = iter
    new_latice[adj_i,adj_j]   = iter
    return new_latice, iter

def getSite(L):
    site_i = rn.randint(0,L-1)
    site_j = rn.randint(0,L-1)
    adj_i  = site_i
    adj_j  = site_j

    side   = rn.randint(0,3)
    if site_i == L-1 and side == 0:
        side = 2
    if site_i == 0 and side == 2:
        side = 0
    if site_j == L-1 and side == 1:
        side = 3
    if site_j == 0 and side == 3:
        side = 1

    if side == 0:
        adj_i += 1
    elif side == 1:
        adj_j += 1
    elif side == 2:
        adj_i -= 1
    elif side == 3:
        adj_j -= 1
    return site_i, site_j, adj_i, adj_j

def annealing(latice,num_dimer,cooling_schedule):
    plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

    FFMpegWriter = ani.writers['ffmpeg']
    metadata = dict(title='annealing', artist='MattDaunt',
                    comment='simulated annealing of dimers')
    writer = FFMpegWriter(fps=20, metadata=metadata)
    outdir = 'outputs/'
    fig = plt.figure()

    length = len(latice)

    ittt = 0
    with writer.saving(fig, outdir + "annealing_dimer_fin1.mp4", 100):
        for temp in cooling_schedule:
            site_i, site_j, adj_i, adj_j = getSite(length)
            if latice[site_i,site_j] == 0 and latice[adj_i,adj_j] == 0:
                num_dimer += 1
                latice[site_i,site_j] = num_dimer
                latice[adj_i,adj_j]   = num_dimer
            elif latice[site_i,site_j] == latice[adj_i,adj_j]:
                dummy = 0
            else:
                new_latice, new_dim = getNewLatice(latice,site_i,site_j,adj_i,adj_j,num_dimer)
                energyDiff = calcEnergy(new_latice) - calcEnergy(latice)
                # print(new_latice==latice)
                if rn.uniform(0,1) < np.exp(-energyDiff/temp):
                    latice    = new_latice
                    num_dimer = new_dim
            if (ittt % 100) == 0:
                plt.imshow(latice)
                writer.grab_frame()
                plt.clf()
            ittt += 1
    return latice

if __name__ == '__main__':
    L = 50
    latice    = np.zeros((L,L))
    num_dimer = 0

    time = np.linspace(0,1,30000)
    Temp = 1
    tau  = .1
    temp_sch   = Temp * np.exp(-time/tau)
    latice_fin = annealing(latice,num_dimer,temp_sch)
    plt.imshow(latice_fin)
    print(calcEnergy(latice_fin))
    # plt.show()
