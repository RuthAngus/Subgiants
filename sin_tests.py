import numpy as np
import matplotlib.pyplot as plt

def fit_sine(x, y, w_true):
    M = np.ones((len(x), 2*len(w_true)+1))
    for i in range(len(w_true)):
        M[:, 2*i] = np.sin(w_true[i]*x)
        M[:, 2*i+1] = np.cos(w_true[i]*x)
    A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
    ys = np.zeros_like(x)
    for i in range(len(w_true)):
        ys += A[2*i]*np.sin(w_true[i]*x) + A[2*i+1]*np.cos(w_true[i]*x)
    ys += A[-1]
    return ys

if __name__ == "__main__":
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

    ys = fit_sine(x, y, w_true)

    plt.clf()
    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.plot(x, ys)
    plt.show()
