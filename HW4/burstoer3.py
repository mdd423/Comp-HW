import numpy as np
import matplotlib.pyplot as plt
import math

# def func(t, y):
    # return -y, y

def Brusslator(t, y):
    a = 1.
    b = 3.
    # print(y)
    dxdt = 1 - ((b+1) * y[0]) + (a * y[0]**2 * y[1])
    dydt =      (b    * y[0]) - (a * y[0]**2 * y[1])
    print('derives: ', dxdt, dydt)
    return np.array([dxdt, dydt])

def modMidSingle(func, time, h, current_val, n):
    #Mod Mid point
    y1 = current_val + (0.5 * h * func(time,current_val))
    y2 = current_val + (h * func(time, y1))
    for iii in range(n-1):
        y1 += h*func(time,y2)
        y2 += h*func(time,y1)
    # if math.isnan(y1[0]) or math.isnan(y2[0]):
    #     print(time, h)
    return 0.5 * (y1 + y2 + 0.5 * h * func(time,y2))

def BerStoerSingle(func, time, current_val, H, delta_t, y_out):
    n = 1
    R1 = np.empty([1,2],float)
    R1[0] = modMidSingle(func,time,H,current_val,n)

    error = 2 * H * delta_t

    while error > H * delta_t:
        if (n > 7):
            # print('entering recursion...', n, '\n')
            y1 = BerStoerSingle(func, time      , current_val, H/2, delta_t, y_out)
            y2 = BerStoerSingle(func, time + H/2, y1         , H/2, delta_t, y_out)
            return y2
        n += 1
        h = H/n

        R2 = R1
        R1 = np.empty([n,2], float)
        R1[0] = modMidSingle(func,time,h,current_val,n)

        for jjj in range(1,n):
            epsilon = (R1[jjj-1] - R2[jjj-1])/((n/(n-1))**(2*jjj)-1)
            R1[jjj] =  R1[jjj-1] + epsilon
        error = abs(epsilon[0])

    y_out[0].append(current_val[0])
    y_out[1].append(current_val[1])
    y_out[2].append(time)

    # current_val = R1[n-1]
    return R1[n-1]

def BerStoer(func, t_i, t_f, N, ivp, delta_t):
    H = (t_f - t_i)/N
    tspace = np.arange(t_i,t_f,H)
    # print(tspace)
    current_val = ivp
    max_n = 8

    y1_out  = [ivp[0]]
    y2_out  = [ivp[1]]
    t_array = [t_i]

    for time in tspace:
        n = 1
        y1 = current_val + 0.5 * H * func(time, current_val)
        y2 = current_val + H * func(time, y1)
        R1 = np.empty([1,2],float)
        R1[0] = 0.5 * (y1 + y2 + 0.5 * H * func(time,y2))

        error = 2*H*delta_t

        max_bool = 0
        while error > H*delta_t:
            n += 1
            y1, y2 = modMidSingle(func, time, H, n, current_val)
            R2 = R1
            R1 = np.empty([n,2], float)
            R1[0] = 0.5 * (y1 + y2 + 0.5 * h * func(time,y2))
            for jjj in range(1,n):
                epsilon = (R1[jjj-1]  - R2[jjj-1])/((n/(n-1))**(2*jjj)-1)
                R1[jjj] =  R1[jjj-1] + epsilon
            error = abs(epsilon[0])
            if (n > max_n):
                max_bool = 1
                print("calling burstoer again...\n")
                y1, y2, ts = BerStoer(func, time, time + H, 2*N, current_val, delta_t)
                y1_out.extend(y1.tolist())
                y2_out.extend(y2.tolist())
                t_array.extend(ts.tolist())
                # current_val[0] = y1
                # current_val[1] = y2
                # t_out = ts
                break
            t_out = time
        current_val = R1[n-1]
        if (max_bool == 0):
            t_array.append(t_out)
            y1_out.append(current_val[0])
            y2_out.append(current_val[1])
    return np.array(y1_out), np.array(y2_out), np.array(t_array)

if __name__ == "__main__":

    delta = 1e-10
    ivp = np.array([0.,0.],float)
    t_i = 0.
    H = 20.

    y_out = [[ivp[0]], [ivp[1]], [t_i]]

    R = BerStoerSingle(Brusslator, t_i, ivp, H, delta, y_out)
    print(y_out)
    # print(y1)
    # print(ts)
    # plt.plot(y_out[1],y_out[2],'o')
    # plt.plot(y_out[0],y_out[2],'d')
    # plt.savefig('Chemical_Reaction')
    # plt.xlabel('Time (unit)')
    # plt.ylabel('Concentration')

    # print(ys)
    # plt.axis('equal')
    # plt.show()
