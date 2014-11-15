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

def fit_sine_err(x, y, yerr, w_true):
    M = np.ones((len(x), 2*len(w_true)+1))
    S = yerr**2 * np.eye(len(yerr))
    for i in range(len(w_true)):
        M[:, 2*i] = np.sin(w_true[i]*x)
        M[:, 2*i+1] = np.cos(w_true[i]*x)
    A = np.linalg.solve(M.T.dot(np.linalg.inv(S)).dot(M),
                        (M.T.dot(np.linalg.inv(S)).dot(y)))
    ys = np.zeros_like(x)
    for i in range(len(w_true)):
        ys += A[2*i]*np.sin(w_true[i]*x) + A[2*i+1]*np.cos(w_true[i]*x)
    ys += A[-1]
    return ys, A

def show_sine(xs, w_true, A):
    ys = np.zeros_like(xs)
    for i in range(len(w_true)):
        ys += A[2*i]*np.sin(w_true[i]*xs) + A[2*i+1]*np.cos(w_true[i]*xs)
    ys += A[-1]
    return ys

if __name__ == "__main__":

    # generate some data
    A_true = [10.1, .4, .6]
    w_true = [20., 10., 5.]
    phi_true = [np.pi/3, np.pi/2., np.pi/6.]
    x = np.linspace(0, 20, 1000)
    y = np.zeros_like(x)
    for i in range(len(A_true)):
        y += A_true[i]*np.sin(w_true[i]*x + phi_true[i])

    plt.clf()
    plt.plot(x, y, 'm')

    err = 1.
    yerr = err * np.random.randn(len(y))
    y += yerr * np.random.randn(len(y))
    ys = fit_sine(x, y, w_true)
    ys2 = fit_sine_err(x, y, yerr, w_true)

    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.plot(x, ys)
    plt.plot(x, ys2, 'g')
    plt.show()
