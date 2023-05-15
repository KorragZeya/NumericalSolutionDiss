import math
import numpy as np
import scipy as sp
import openpyxl

book = openpyxl.Workbook()
sheet = book.active


def q(k, dt):
    if (k*dt//43200) % 2 == 0:
        return 1500*np.sin((2 * 3.14 / 86400) * k * dt)
    else:
        return 0


def numerical_solution():
    L = 3
    Nx = 100
    end_index = Nx - 1
    t = 2.5 * 86000
    dt = 100
    Nt = int(t / dt)
    dx = L / Nx
    Initial = np.full(Nx, 288)
    lambdas_h = np.full(Nx + 1, 0.7)
    Cs = 2006200
    Cw = 4200000
    F = 0

    Tair = lambda k: 288 + 5 * np.sin((2 * 3.14 / 86400) * k * dt)

    M = np.zeros((Nx, Nx), dtype=float)

    alphar = 0.6

    for k in range(Nt):
        for i in range(Nx):
            a = -(dt / dx ** 2) * lambdas_h[i]
            c = -(dt / dx ** 2) * lambdas_h[i + 1]
            b = Cs - c - a
            if i == 0:
                M[0][0] = b - a*(2*dx/lambdas_h[0])*0.9*5.7*10**(-8)*Initial[0]**3
                M[0][1] = a + c
            elif i == Nx - 1:
                M[end_index][end_index - 1] = a
                M[end_index][end_index] = b + c
            else:
                M[i][i - 1] = a
                M[i][i] = b
                M[i][i + 1] = c

        D = np.array([Cs * Initial[0] - a*(2*dx/lambdas_h[0])*0.8*0.9*5.7*10**(-8)*Tair(k)**4 - a*(2*dx/lambdas_h[0])*q(k, dt) if i == 0
                      else Cs * Initial[i] for i in range(Nx)], dtype=float)

        Initial = np.copy(sp.sparse.csr_matrix(sp.linalg.inv(M)).dot(D))

        print(Initial[0]-273)

    for i in range(Nx):
        sheet.cell(row=i + 1, column=1).value = float(Initial[i]) - 273
        sheet.cell(row=i + 1, column=2).value = L - i * dx

    book.save("Numerical_solution.xlsx")
    book.close()


numerical_solution()
