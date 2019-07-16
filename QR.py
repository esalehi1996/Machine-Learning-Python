import numpy as np
np.set_printoptions(precision=9)


def QR(M,n):
    u = np.copy(M).transpose()[:][0]
    T = np.identity(n, dtype=float)
    for i in range(n - 1):
        u[0] = u[0] - np.linalg.norm(u)
        if (np.linalg.norm(u)==0):
            u = np.zeros(shape=(1, len(u)))
        else:
            u = u / np.linalg.norm(u)
        A = np.identity(n - i, dtype=float) - 2 * np.matmul(u.reshape(n - i, 1), u.reshape(1, n - i))
        B = np.zeros((n, n), dtype=float)
        for j in range(i):
            B[j, j] = 1
        B[i:, i:] = A
        T = np.matmul(B, T)
        u = np.matmul(T, M)[i + 1:, i + 1]
    return T.transpose(),np.matmul(T, M)

n = int(input())

M = np.zeros((n, n), dtype=float)

for i in range(n):
    M[i] = [float(i) for i in input().split(' ')]


U = np.identity(n , dtype=float)
S = M.transpose()
V = np.identity(n , dtype=float)

loopindex = 0
tol = 1e-10
Error = 1
while ((Error > tol) & (loopindex < 100*n)):
    Q , S = QR(S.transpose(),n)
    U = np.matmul(U,Q)
    Q , S = QR(S.transpose(),n)
    V = np.matmul(V,Q)
    
    F = 0
    for i in range(n):
        F = F + S[i,i]*S[i,i]
    E = np.linalg.norm(S)*np.linalg.norm(S) - F
    F = np.sqrt(F)
    E = np.sqrt(E)
    if (F == 0):
        F=1
    Error = E/F
    loopindex = loopindex + 1
    
for i in range(n):
    if (S[i,i]<0):
        U[:,i] = -U[:,i]
    S[i,i] = np.abs(S[i,i])

for i in U:
    for j in i:
        print('%.9f' % j, end=' ')
    print()

for i in S:
    for j in i:
        print('%.9f' % j, end=' ')
    print()

for i in V:
    for j in i:
        print('%.9f' % j, end=' ')
    print()















