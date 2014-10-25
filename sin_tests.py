import numpy as np
import matplotlib.pyplot as plt

# generate some data
A_true = [.1, .4, .6]
w_true = [2., 1., .5]
phi_true = [np.pi/3, np.pi/2., np.pi/6.]

x = np.linspace(0, 20, 1000)
y = np.zeros_like(x)
for i in range(len(A_true)):
    y += A_true[i]*np.sin(w_true[i]*x + phi_true[i])
yerr = np.ones_like(x)*0.01
y += 0.1 * np.random.randn(len(y))

M = np.ones((2*len(w_true)+1, len(x)))
for i in range(len(w_true)):
    M[2*i, :] = np.sin(w_true[i]*x)
    M[2*i+1, :] = np.cos(w_true[i]*x)
M = M.T

A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))

print np.shape(A)
ys = np.zeros_like(x)
for i in range(len(w_true)):
    ys += A[i*i]*np.sin(w_true[i]*x) + A[2*i+1]*np.cos(w_true[i]*x)

plt.clf()
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(x, ys)
plt.show()
