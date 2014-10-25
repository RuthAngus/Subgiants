import numpy as np
import matplotlib.pyplot as plt

# generate some data
# A_true = [.1, .4, .6]
# w_true = [2., 1., .5]
# phi_true = [np.pi/3, np.pi/2., np.pi/6.]

A_true = [.1]
w_true = [2.]
phi_true = [np.pi/3.]
x = np.linspace(0, 20, 1000)
y = np.zeros_like(x)
for i in range(len(A_true)):
    y += A_true[i]*np.sin(w_true[i]*x + phi_true[i])
yerr = np.ones_like(x)*0.01
y += 0.1 * np.random.randn(len(y))

M = np.vstack((np.sin(w_true*x), np.cos(w_true*x))).T
A = np.linalg.inv(np.dot(M.T, M)) * np.dot(M.T, y)

print np.shape(A)
ys = np.zeros_like(x)
for i in range(len(A)):
    ys += A[i][0]*np.sin(w_true[0]*x) + A[i][1]*np.cos(w_true[0]*x)

plt.clf()
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
plt.plot(x, ys)
plt.show()

# M = np.vstack((np.sin(w_true[0]*x), np.cos(w_true[0]*x),
#                np.sin(w_true[1]*x), np.cos(w_true[1]*x),
#                np.sin(w_true[2]*x), np.cos(w_true[2]*x))).T
# M = np.zeros((len(w_true), 2*len(w_true))
# for i in range(len(w_true)):
#     M[i][ =
