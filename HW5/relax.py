import numpy as np
import matplotlib.pyplot as plt

t_i = 0
t_f = 10
n = 100
h = (t_f - t_i)/n

g = 9.8
psi = np.zeros((n,2))
stop = 1000

dotheloop = 1
iter = 0
indiv_error = 1e-10
tot_error = n * indiv_error
step = 0
while dotheloop:
    psi[1:n-1,step] = 1/2 * (psi[0:n-2,(step-1)%2] + psi[2:n,(step-1)%2] + (g*h**2))
    iter += 1
    step += 1
    step = step % 2
    sum_psi = np.sum(psi,axis=0)
    difference = abs(sum_psi[1] - sum_psi[0])
    if difference < tot_error:
        print(iter)
        break

plt.plot(psi)
plt.xlabel('Time (s)')
plt.ylabel('Y distance (m)')
plt.savefig('falling_obj.png')
