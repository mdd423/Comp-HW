import numpy as np
import matplotlib.pyplot as plt
import time

def approxEqui(x1, x2, epsilon):
    if (((x1 + epsilon) > x2) and ((x1 - epsilon) < x2)):
        return True
    else:
        return False

def func(t, y):
    return -t

def gravitation(t, y):
    A, M, G, B = 1., 1., 1., 1.
    r_mag = np.sqrt(y[0]**2 + y[1]**2)
    v_mag = np.sqrt(y[2]**2 + y[3]**2)
    v_df = -A / (v_mag**3 + B)
    # v_df = 0
    return np.array([y[2], y[3], (-(G * M )/r_mag**3)*y[0] + (v_df * y[2]), (-(G * M )/r_mag**3)*y[1] + (v_df*y[3])])

def rk45(func, tspace, ivp):

    y_out = np.zeros((tspace.shape[0], ivp.shape[0]))
    y_out[0,:] = ivp
    N = len(tspace)
    delta_t = (tspace[0] - tspace[-1])/N

    for iii in range(1,N):
        k1 = delta_t * func(tspace[iii], y_out[iii-1,:])
        k2 = delta_t * func(tspace[iii] + delta_t/2, y_out[iii-1,:] + k1/2)
        k3 = delta_t * func(tspace[iii] + delta_t/2, y_out[iii-1,:] + k2/2)
        k4 = delta_t * func(tspace[iii] + delta_t, y_out[iii-1,:] + k3)

        y_out[iii] = y_out[iii-1] + (1/6) * (k1 + 2*k2 + 2*k3 + k4)

    return y_out

def rk_1step(func,ivp,time,dt):
    k1 = dt * func(time, ivp)
    k2 = dt * func(time + dt/2, ivp + k1/2)
    k3 = dt * func(time + dt/2, ivp + k2/2)
    k4 = dt * func(time + dt, ivp + k3)
    return ivp + (1/6) * (k1 + 2*k2 + 2*k3 + k4)

def rk45varystep(func, t_i, t_f, ivp):
    allowed_error = 1e-4
    time_current = t_i
    time_array = np.array([t_i])
    delta_t = (t_f - t_i)/2
    y_out = np.zeros([1,ivp.shape[0]])
    y_out[0,:] = ivp
    iter = 0
    cnt = 0
    max_iter = 10000000

    start_time = time.time()
    current_ten = 0

    grav_file = open('grav_file2.txt','w')
    while (time_current < t_f):
        elapsed_time = time.time() - start_time
        if (elapsed_time//600) > current_ten:
            print('passing ',elapsed_time//600 * 10, ' min')
        current_ten = elapsed_time//600

        cnt += 1
        step_1 = rk_1step(func, rk_1step(func, y_out[iter,:], time_current, delta_t), time_current, delta_t)
        step_2 = rk_1step(func, y_out[iter,:], time_current, 2*delta_t)
        # trying to side step cancellation errors
        if (approxEqui(np.sum(step_1), np.sum(step_2), 1e-6)):
            diff_calc = 1e-6
        else:
            diff_calc = np.sum(np.abs(step_1 - step_2))

        rho = (30 * delta_t * allowed_error)/diff_calc
        if (rho >= 1):
            grav_file.write(str(time_current) + '\t' + str(step_1))
            y_out = np.concatenate((y_out,[step_1]))
            # y_out.append(step_1)
            time_current += 2 * delta_t
            # time_array.append(time_current)
            time_array = np.append(time_array,time_current)
            if (rho**(1/4) > 2):
                delta_t *= 2
            else:
                delta_t *= rho**(1/4)
            iter += 1
        else:
            delta_t *= rho**(1/4)
        # if (cnt > max_iter):
        #     print("count exceeded")
        #     break
    grav_file.close()
    return y_out, time_array


if __name__ == "__main__":
    ivp = np.array([1,0,0,1])
    t_i = 0
    t_f = 15

    time_unit_conv = (4.300e-3)*(100**(-3))*(3.2408e-14)**2 / (1e-8 * (3.171e-8)*2)
    # ys, ts = rk45varystep(func, t_i, t_f, 10)
    ys, ts = rk45varystep(gravitation, t_i, t_f, ivp)
    plt.plot(ys[:,0],ys[:,1])
    plt.title('Orbit')
    plt.axis('equal')
    plt.savefig('grav_orbit_3_vdf.png')
    plt.show()

    plt.clf()

    plt.plot(np.sqrt(1/time_unit_conv) * ts, (ys[:,0])**2 + (ys[:,1])**2)
    plt.ylabel('radius (100 Mpc)')
    plt.xlabel('time (yr)')
    plt.savefig('grav_rad_3_vdf.png')
    plt.show()
