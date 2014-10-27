import numpy as np
import matplotlib.pyplot as plt

def fit_sine(x, y, w):
    m = np.ones((len(x), 2*len(w)+1))
    for i in range(len(w)):
        m[:, 2*i] = np.sin(w[i]*x)
        m[:, 2*i+1] = np.cos(w[i]*x)
    a = np.linalg.solve(np.dot(m.T, m), np.dot(m.T, y))
    ys = np.zeros_like(x)
    for i in range(len(w)):
        ys += a[2*i]*np.sin(w[i]*x) + a[2*i+1]*np.cos(w[i]*x)
    ys += a[-1]
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
