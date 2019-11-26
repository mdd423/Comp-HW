import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sparse

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

    ones = np.ones(M)
    d2_vec = [ones,-2*ones,ones]
    D2 = sparse.spdiags(d2_vec,[-1,0,1],M,M)
    I  = sparse.spdiags(ones,0,M,M)
    Lap = sparse.kron(D2,I) + sparse.kron(I,D2)

    stop = 1000

    dotheloop = 1
    iter = 0
    indiv_error = 1e-10
    tot_error = M**2 * indiv_error

    charge_vec = np.reshape(charge_dist, (M*M))

    step = 0
    phi = np.zeros((M*M,2))
    weight = 0.1
    for hhh in range(stop):
        for iii in range(M*M):
            print("element: ", iii)
            temp_lower = 0
            temp_upper = 0
            for ggg in range(M*M):
                if ggg < iii:
                    temp_lower += Lap[iii,ggg] * phi[ggg,step]
                elif ggg > iii:
                    temp_upper += Lap[iii,ggg] * phi[ggg,(step-1)%2]
            phi[iii,step] = (1 - weight) * phi[iii, (step-1)%2] + weight/Lap[iii,iii] * (charge_vec[iii] - temp_lower - temp_upper)
        print("iter: ", hhh)
        step += 1
        step %= 2

    phi_final = np.reshape(phi[:,:,step],(M,M))
    plt.imshow(phi_final, cmap='gray')
    plt.savefig('potential_plot_2.png')
    plt.show()
    plt.clf()
