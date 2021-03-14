import numpy as np

def f(x,y):
    return -1 * (np.pi ** 2) * (x**2 + y**2) * np.exp((np.sin(np.pi * x * y)) ** 2) * (2 * np.cos(2 * np.pi * x * y) + (np.sin(2 * np.pi * x * y))**2)

def u(x,y):
    return np.exp((np.sin(np.pi * x * y)) ** 2)

def m1(y):
    return 1

def m2(y):
    return np.exp((np.sin(np.pi * y)) ** 2)

def m3(x):
    return 1

def m4(x):
    return np.exp((np.sin(np.pi * x)) ** 2)

def trueMatrix(n,m):
    U = np.zeros((m + 1, n + 1))
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b - a) / n
    k = (d - c) / m
    for j in range(m + 1):
        for i in range(n + 1):
            U[i][j] = u(a + i * h, c + j * k)
    return U

def initialMatrixF(n, m):
    F = np.zeros((m + 1, n + 1))
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b - a) / n
    k = (d - c) / m
    for j in range(1,m):
        for i in range(1,n):
            F[i][j] = f(a + i * h, c + j * k)
#    print(F)
#    print("up F")
    return F

def initialMatrixV(n, m):
    V = np.zeros((m + 1, n + 1))
    print(V)
    print(n,m)
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b-a)/n
    k = (d-c)/m
    for j in range(m + 1):
        V[0][j] = m1(c + j * k)
    for j in range(m + 1):
        V[n][j] = m2(c + j * k)
    for i in range(n + 1):
        V[i][0] = m3(a + i * h)
    for i in range(n + 1):
        V[i][m] = m4(a + i * h)
    return V

def residual(n, m, V, F):
    res = np.zeros((m + 1, n + 1))
    h_quad = -1 * n * n
    k_quad = -1 * m * m
    A = -2 * (h_quad + k_quad)
    for j in range(1,m):
        for i in range(1,n):
            res[i][j] = A * V[i][j] + h_quad * (V[i - 1][j] + V[i + 1][j])
            res[i][j] = res[i][j] + k_quad * (V[i][j + 1] + V[i][j - 1]) - F[i][j]
    return res

def tsCounter(n, m, V, F, re):
    scal1, scal2 = 0, 0
    h_quad = -1 * n * n
    k_quad = -1 * m * m
    A = -2 * (h_quad + k_quad)
    for j in range(1, m):
        for i in range(1,n):
            up, down, left, right = 1, 1, 1, 1
            if (i == 1):
                left = 0
            if (i == n - 1):
                right = 0
            if (j == 1):
                down = 0
            if (j == m - 1):
                up = 0
            temp = A * re[i][j] + h_quad * (left * re[i - 1][j] + right * re[i + 1][j])
            temp += k_quad * (up * re[i][j + 1] + down * re[i][j - 1])
            scal1 += temp * re[i][j]
            scal2 += temp * temp
    return scal1 / scal2

def methodMR(n, m, SMax, eps, S):
    eps_max = 0
    eps_cur = 0
    S[0] = 0
    V = initialMatrixV(n,m)
    F = initialMatrixF(n,m)
    v_old = 0
    v_new = 0
    flag = False
    while(flag != True):
        eps_max = 0
        re = residual(n,m,V,F)
        ts = tsCounter(n,m,V,F,re)
#        print("re")
#        print(re)
#        print(ts)
#        print("V")
#        print(V)
        temp = 0
        for j in range(1, m):
            for i in range(1, n):
                v_old = V[i][j]
                v_new = V[i][j] - ts * re[i][j]
                temp += re[i][j] * re[i][j]
                eps_cur = np.fabs(v_old - v_new)
                if (eps_cur > eps_max):
                    eps_max = eps_cur
                V[i][j] = v_new
        S[0] = S[0] + 1
        S[2] = np.sqrt(temp)

        if ((eps_max < eps) or (S[0] >= SMax)):
            flag = True
    S[1] = eps_max
    return V

n = 4
m = 3
F = np.zeros((m + 1, n + 1))
F[1][1] = 1
F[0][1] = -1
F[1][0] = -2
print(F)