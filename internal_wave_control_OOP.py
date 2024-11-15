# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 20:08:39 2024

@author: dangt
"""

import numpy as np
from scipy.sparse import spdiags, eye
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt


class CoupledWaveFVM:
    def __init__(self, a, b, T0, T, n, m, alpha, gamma, mu, bl, br):
        self.a = a
        self.b = b
        self.T0 = T0
        self.T = T
        self.n = n
        self.m = m
        self.alpha = alpha
        self.gamma = gamma
        self.mu = mu
        self.bl = bl
        self.br = br
        self.dx = (b - a) / (n - 1)
        self.dt = (T - T0) / (m - 1)
        self.x = np.linspace(a, b, n)
        self.t = np.linspace(T0, T, m)

    def solve(self, y0, y0t, z, yhat, lam_T):
        s = self.dt ** 2 / self.dx ** 2
        s1 = self.dt / self.dx

        A = spdiags([s * np.ones(self.n - 2), 2 * (1 - s) * np.ones(self.n - 2), s * np.ones(self.n - 2)],
                    [-1, 0, 1], self.n - 2, self.n - 2)
        C = spdiags([s1 / 2 * np.ones(self.n - 2), -s1 * np.ones(self.n - 2), s1 / 2 * np.ones(self.n - 2)],
                    [-1, 0, 1], self.n - 2, self.n - 2)

        Yhat = np.zeros((2 * (self.n - 2) * self.m, 1))
        Yhat[(self.n - 2):2 * (self.n - 2), 0] = self.dt * y0t[1:-1]
        Yhat[:self.n - 2, 0] = y0[1:-1]

        B = eye(2 * (self.n - 2) * self.m)
        Z = np.zeros((2 * (self.n - 2) * self.m, 1))
        Z[((2 * self.m - 1) * (self.n - 2)):, 0] = lam_T[1:-1]
        Z[((2 * self.m - 2) * (self.n - 2)):((2 * self.m - 1) * (self.n - 2)), 0] = \
            self.dt * z[1:-1] + 0.5 * self.dt ** 2 * yhat[1:-1, -1]

        for i in range(self.m):
            k1 = i + 1
            k2 = i + 1 + self.m
            if k1 > 2:
                B[(k1 - 1) * (self.n - 2):(k1 * (self.n - 2)),
                  (k1 - 2) * (self.n - 2):(k1 - 1) * (self.n - 2)] = eye(self.n - 2) + self.dt / self.dx * C
            if k2 < (2 * self.m - 1):
                B[(k2 - 1) * (self.n - 2):(k2 * (self.n - 2)),
                  (k2 - self.m - 1) * (self.n - 2):(k2 - self.m + 1) * (self.n - 2)] = -self.dt ** 2 * eye(self.n - 2)
                B[(k2 - 1) * (self.n - 2):(k2 * (self.n - 2)),
                  k2 * (self.n - 2):(k2 + 1) * (self.n - 2)] = eye(self.n - 2) + self.dt / self.dx * C
            if k2 == (2 * self.m - 1):
                B[(k2 - 1) * (self.n - 2):(k2 * (self.n - 2)),
                  (k2 - self.m - 1) * (self.n - 2):(k2 - self.m + 1) * (self.n - 2)] = \
                    -0.5 * self.dt ** 2 * eye(self.n - 2) - self.dt * self.gamma * eye(self.n - 2)
                B[(k2 - 1) * (self.n - 2):(k2 * (self.n - 2)),
                  k2 * (self.n - 2):(k2 + 1) * (self.n - 2)] = eye(self.n - 2) + self.dt / self.dx * C

        M = eye(2 * (self.n - 2) * self.m) - B
        GrandVector = spsolve(M, Yhat + Z)

        Y = np.zeros((self.n, self.m))
        Lambda = np.zeros((self.n, self.m))
        for i in range(self.m):
            Y[1:-1, i] = GrandVector[i * (self.n - 2):(i + 1) * (self.n - 2), 0]
            Lambda[1:-1, i] = GrandVector[((i + self.m) * (self.n - 2)):((i + self.m + 1) * (self.n - 2)), 0]

        return Y, Lambda


class ErrorAnalysis:
    @staticmethod
    def compute_error(Y, y_exact, Lambda, L_exact, dt, dx):
        error_y = (dt * dx) ** (1 / 2) * np.linalg.norm(Y - y_exact)
        error_max = np.max(np.abs(Y - y_exact))
        error_control = (dt * dx) ** (1 / 2) * np.linalg.norm(Lambda - L_exact)
        return error_y, error_max, error_control

a = 0
b = 1
T0 = 0
TT = np.arange(1, 51, 1)
alpha = 1
gamma = 0
mu = 1
error = []
error_control = []
error_max = []
relative = []
N = [51, 101, 201, 401]
M = [51, 101, 201, 401]

for l in range(len(N)):
    T = 1
    n = N[l]
    m = M[l]
    dx = (b - a) / (n - 1)
    dt = T / (m - 1)
    bl = np.zeros(m)
    br = np.zeros(m)
    x = np.linspace(a, b, n)
    t = np.linspace(0, T, m)

    # Define the y hat solution
    yhat = np.zeros((n, m))
    z = np.zeros(n)  # gamma*y(x,T) - lambda_t(x,T)

    # Test for example y(x,t) = sin(pi*x)*(t-T)
    y0 = -np.sin(np.pi * x) * T
    y0t = np.sin(np.pi * x)

    for i in range(len(x)):
        # Test for example y(x,t) = sin(pi*x)*(t-T)
        yhat[i, :] = (alpha * np.pi**4 + 1) * np.sin(np.pi * x[i]) * (t - T)

    # Test for example y(x,t) = sin(pi*x)*(t-T)
    z = -alpha * np.pi**2 * np.sin(np.pi * x)

    # Define the exact solution y
    yexact = np.zeros((n, m))
    Lexact = np.zeros((n, m))
    for i in range(len(x)):
        # Exact solution for example sin(pi*x)*(t-T)
        yexact[i, :] = np.sin(np.pi * x[i]) * (t - T)
        Lexact[i, :] = alpha * np.pi**2 * np.sin(np.pi * x[i]) * (t - T)

    lam_T = Lexact[:, -1]

    # Instantiate the CoupledWaveFVM class and solve
    solver = CoupledWaveFVM(a, b, T0, T, n, m, alpha, gamma, mu, bl, br)
    Y, Lambda = solver.solve(y0, y0t, z, yhat, lam_T)

    # Visualization
    for j in range(len(t)):
        plt.subplot(2, 1, 1)
        plt.plot(x, Y[:, j])
        plt.title('Y', fontsize=30)
        plt.axis([0, 1, -1, 1])
        plt.subplot(2, 1, 2)
        plt.plot(x, Lambda[:, j])
        plt.title('\Lambda', fontsize=30)
        plt.pause(0.5)
    plt.show()

    # Error Analysis
    err, err_max, err_control = ErrorAnalysis.compute_error(Y, yexact, Lambda, Lexact, dt, dx)
    error.append(err)
    error_max.append(err_max)
    error_control.append(err_control)

# Plot the order of convergence
plt.figure()
plt.loglog(N, np.power(N, -2), N, error, '-x', N, error_max, '-+', N, error_control, '-o', linewidth=2.0)
plt.legend(['2x', '||y_d - y_{exact}||_{l^{2}}', '||y_d - y_{exact}||_{l^{\infty}}', '||\lambda_d - \lambda_{exact}||_{l^{2}}'], fontsize=30)
plt.xlabel('N', fontsize=30)
plt.ylabel('error', fontsize=30)
plt.gca().set_fontsize(30)
plt.show()

