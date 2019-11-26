import numpy as np
import matplotlib.pyplot as plt

def readIn(file):
    f = open(file,"r")
    lines = f.readlines()
    result_x = []
    result_y = []
    iter = 0
    for x in lines:
        temp = x.split()
        result_x.append(float(temp[0]))
        result_y.append(float(temp[1]))
        # if (iter % 100) == 0:
        #     print(float(x.split(' ')))
        iter += 1
    f.close()
    new_list = [result_x, result_y]
    # print(new_list)
    return np.asarray(new_list, dtype=float)

def cloudincell(particles, M):
    size = len(particles)

    charge_array = np.zeros((M,M))
    charge_dist  = np.zeros((M,M))

    for iii in range(size):
        charge_array[int(particles[iii,0]), int(particles[iii,1])] += 1
        if particles[iii,0] < 0.5:
            if particles[iii,1] < 0.5:
                 a = particles[iii,0] % 1
                 b = particles[iii,1] % 1
                 charge_dist[int(particles[iii,0]), int(particles[iii,1])]     += (0.5 + a) * (0.5 + b)
                 charge_dist[int(particles[iii,0]) - 1, int(particles[iii,1])] += (0.5 - a) * (0.5 + b)
                 charge_dist[int(particles[iii,0]), int(particles[iii,1]) - 1] += (0.5 + a) * (0.5 - b)
                 charge_dist[int(particles[iii,0]) - 1, int(particles[iii,1]) - 1] += (0.5 - a) * (0.5 - b)
            else:
                a = particles[iii,0] % 1
                b = 1 - (particles[iii,1] % 1)
                charge_dist[int(particles[iii,0]), int(particles[iii,1])]     += (0.5 + a) * (0.5 + b)
                charge_dist[int(particles[iii,0]) - 1, int(particles[iii,1])] += (0.5 - a) * (0.5 + b)
                charge_dist[int(particles[iii,0]), int(particles[iii,1]) + 1] += (0.5 + a) * (0.5 - b)
                charge_dist[int(particles[iii,0]) - 1, int(particles[iii,1]) + 1] += (0.5 - a) * (0.5 - b)

        else:
            if particles[iii,1] < 0.5:
                a = 1 - (particles[iii,0] % 1)
                b = particles[iii,1] % 1
                charge_dist[int(particles[iii,0]), int(particles[iii,1])]     += (0.5 + a) * (0.5 + b)
                charge_dist[int(particles[iii,0]) + 1, int(particles[iii,1])] += (0.5 - a) * (0.5 + b)
                charge_dist[int(particles[iii,0]), int(particles[iii,1]) - 1] += (0.5 + a) * (0.5 - b)
                charge_dist[int(particles[iii,0]) + 1, int(particles[iii,1]) - 1] += (0.5 - a) * (0.5 - b)
            else:
                a = 1 - (particles[iii,0] % 1)
                b = 1 - (particles[iii,1] % 1)
                charge_dist[int(particles[iii,0]), int(particles[iii,1])]     += (0.5 + a) * (0.5 + b)
                charge_dist[int(particles[iii,0]) + 1, int(particles[iii,1])] += (0.5 - a) * (0.5 + b)
                charge_dist[int(particles[iii,0]), int(particles[iii,1]) + 1] += (0.5 + a) * (0.5 - b)
                charge_dist[int(particles[iii,0]) + 1, int(particles[iii,1]) + 1] += (0.5 - a) * (0.5 - b)
    return charge_dist

def gaussSidel(weight, charge_dist):
    phi  = np.zeros((M+2,M+2,2))

    dotheloop = 1
    iter = 0
    indiv_error = 1e-10
    tot_error = M**2 * indiv_error
    step = 0

    while dotheloop:
        for iii in range(1,M+1):
            for jjj in range(1,M+1):
                phi[iii, jjj, (step+1)%2] = (-weight * phi[iii,jjj,step]) + ((1+weight)/4. * (phi[iii-1, jjj, (step+1)%2] + phi[iii+1, jjj, step] + phi[iii, jjj-1, (step+1)%2] + phi[iii, jjj+1, step] + charge_dist[iii-1, jjj-1]))
        iter += 1
        step += 1
        step %= 2
        sum_phis = np.sum(np.sum(phi,axis=0),axis=0)

        difference = abs(sum_phis[1] - sum_phis[0])

        if (iter % 100) == 0:
            print(iter)
            print(difference)
        if difference < tot_error:
            return iter

def goldenRaio(gauss, a, b, charge_dist):
    gr = (np.sqrt(5) + 1) / 2

    ab_array = [[a,b]]

    error = 1e-3

    h = (b - a) / gr
    c = b - h
    d = a + h
    while abs(c - d) > error:
        if gauss(c, charge_dist) < gauss(d, charge_dist):
            b = d
        else:
            a = c

        h = (b - a) / gr
        c = b - h
        d = a + h
        ab_array.append([a,b])
    return np.asarray(ab_array)

if __name__ == "__main__":
    particles = np.transpose(readIn('particles.dat'))
    M = 100

    charge_dist = cloudincell(particles, M)
    charge_dist *= (1.602e-19/8.85e-12)

    ab_array = goldenRaio(gaussSidel, .95, .98, charge_dist)
    print((ab_array[-1,0] + ab_array[-1,1]) / 2)

    # plt.plot(ab_array[:,0])
    # plt.plot(ab_array[:,1])
    # plt.savefig('optimal_weight.png')
    # plt.show()
