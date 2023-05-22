import math
import numpy as np
import scipy as sp
import openpyxl

book = openpyxl.Workbook()
sheet = book.active


"""
это функция радиационного солнечного потока. 
Если время соответствует ночи, то поток равен 0.
"""
def q(k, dt):
    if (k*dt//43200) % 2 == 0:
        qin = 700*np.sin((2 * 3.14 / 86400) * k * dt)
        if qin < 0:
            qin = 0
        return qin
    else:
        return 0



def numerical_solution():
    L1, L2, L3, L4 = 0.10, 0.20, 1.20, 1.00
    L = L1 + L2 + L3 + L4
    # L = 3
    t = 1.4 * 86400
    dt = 25
    Nt = int(t / dt)
    dx = 0.01
    Nx = int(L/dx)
    end_index = Nx - 1
    Initial = np.full(Nx, 273)

    λ1, λ2, λ3, λ4 = 0.7, 0.9, 0.7, 1.7
    # lambdas_h = np.full(Nx + 1, 0.7)
    lambda_i = [0.07] + [λ1 for _ in range(int(L1 / dx))] + \
               [λ2 for _ in range(int(L2 /dx))] + \
               [λ3 for _ in range(int(L3 / dx))] + \
               [λ4 for _ in range(int(L4 / dx))] + [λ4]

    lambdas_h = np.array([2 * lambda_i[i] * lambda_i[i + 1] / (lambda_i[i] + lambda_i[i + 1]) for i in range(Nx + 1)])

    Cs = 2006200
    Cw = 4200000
    F = 0

    Tair = lambda k: 288 + 5 * np.sin((2 * 3.14 / 86400) * k * dt)

    M = np.zeros((Nx, Nx), dtype=float)

    for k in range(Nt):
        alpha = 10 * (Initial[0] - Tair(k))**0.33 if Initial[0] > Tair(k) else 10 * (Tair(k) - Initial[0])**0.33
        for i in range(Nx):
            a = -(dt / dx ** 2) * lambdas_h[i]
            c = -(dt / dx ** 2) * lambdas_h[i + 1]
            b = Cs - c - a
            if i == 0:
                M[0][0] = b + a
                M[0][1] = c
            elif i == Nx - 1:
                M[end_index][end_index - 1] = a
                M[end_index][end_index] = b + c
            else:
                M[i][i - 1] = a
                M[i][i] = b
                M[i][i + 1] = c
        print(Initial[0] - 273)
        D = np.array([Cs * Initial[0]
                      - a * (dx / lambdas_h[0]) * 0.83 * 0.8 * 5.7 * 10 ** (-8) * Tair(k) ** 4
                      + a * (dx / lambdas_h[0]) * 0.83 * 5.7 * 10 ** (-8) * Initial[0] ** 4
                      - a * (dx/lambdas_h[0]) * q(k, dt)
                      + a * (dx / lambdas_h[0]) * alpha * (Initial[0] - Tair(k))
                      if i == 0 else Cs * Initial[i] for i in range(Nx)], dtype=float)

        Initial = np.copy(sp.sparse.csr_matrix(sp.linalg.inv(M)).dot(D))

    print(type(float(Initial[2])))

    for i in range(Nx):
        sheet.cell(row=i + 1, column=1).value = float(Initial[i]) - 273
        sheet.cell(row=i + 1, column=2).value = L - i * dx

    book.save("Numerical_solution.xlsx")
    book.close()


numerical_solution()
