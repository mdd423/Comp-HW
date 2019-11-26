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

# def cloudincell(array, ):


if __name__ == "__main__":
    particles = np.transpose(readIn('particles.dat'))
    size = len(particles)
    # floor_parts = int(np.floor(particles))
    M = 100
    charge_array = np.zeros((M,M))
    charge_dist  = np.zeros((M,M))
    # print(charge_dist.shape)

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

    # plt.imshow(charge_dist, cmap='gray')
    # plt.savefig('charges.png')
    # plt.show()
    # plt.clf()

    # plt.imshow(charge_array, cmap='gray')
    # plt.savefig('charge_array.png')
    # plt.show()
    # plt.clf()

    phi  = np.zeros((M+2,M+2,2))
    stop = 1000
    # test_arr = np.arange(10)
    # print(phi[M,M])
    dotheloop = 1
    iter = 0
    indiv_error = 1e-10
    tot_error = M**2 * indiv_error
    step = 0
    while dotheloop:
        phi[1:(M+1), 1:(M+1), (step+1)%2] = 1./4. * (phi[1:(M+1), 2:, step] + phi[0:M, 1:(M+1), step] + phi[1:(M+1), 0:M, step] + phi[2:, 1:(M+1), step] + charge_dist)
        iter += 1
        step += 1
        step = step % 2
        sum_phis = np.sum(np.sum(phi,axis=0),axis=0)
        # print(sum_phis)
        difference = abs(sum_phis[1] - sum_phis[0])
        if difference < tot_error:
            print(iter)
            break

    plt.imshow(phi[:,:,(step-1) % 2], cmap='gray')
    plt.savefig('potential_plot.png')
    plt.show()
    plt.clf()
