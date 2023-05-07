import math
import numpy as np
import scipy as sp
import openpyxl

book = openpyxl.Workbook()

sheet = book.active

def numerical_solution():
    L1, L2, L3, L4 = 0.10, 0.20, 0.20, 1.50
    L = L1 + L2 + L3 + L4
    λ1, λ2, λ3, λ4 = 0.4, 0.5, 0.9, 0.5
    Nx = 100
    end_index = Nx - 1
    t = 2.3 * 86000
    dt = 25
    Nt = int(t/dt)
    dx = L/Nx
    Initial = np.full(Nx, 288)

    lambda_i = [0.025]+[λ1 for i in range(int(L1/L*Nx))]+\
               [λ2 for i in range(int(L2/L*Nx))]+\
               [λ3 for i in range(int(L3/L*Nx))]+\
               [λ4 for i in range(int(L4/L*Nx))]+[λ4]

    # lambdas_h = np.full(Nx+1, 0.7)
    # lambdas_h = np.array([0.04]+[0.1 for i in range(Nx)])
    lambdas_h = np.array([2*lambda_i[i]*lambda_i[i+1]/(lambda_i[i]+lambda_i[i+1]) for i in range(Nx+1)])

    Cs = 2006200
    Cw = 4200000
    F = 0

    q = lambda k: 700*np.sin((2*3.14/86400)*k*dt)

    Tair = lambda k: 288 + 0*np.sin((2*3.14/86400)*k*dt)

    M = np.zeros((Nx, Nx), dtype=float)


    for k in range(Nt):
        alphac = 10 * math.pow(abs(Initial[0] - Tair(k)), 0.33)
        alphar = 0.9 * 5.67 * 10 ** (-8) * \
                 (Initial[0] ** 2 - math.pow(0.8, 0.5) * Tair(k) ** 2) * \
                 (Initial[0] - math.pow(0.8, 0.25) * Tair(k))
        alpha = alphac + alphar

        D = np.array([Cs*Initial[0] + -(dt/dx**2)*lambdas_h[0]*2*dx*alpha/lambdas_h[0]*Tair(k) - -(dt/dx**2)*lambdas_h[0]*2*dx/lambdas_h[0]*q(k) if i == 0
                      else Cs*Initial[i] for i in range(Nx)], dtype=float)

        for i in range(Nx):
                a = -(dt/dx**2)*lambdas_h[i]
                c = -(dt/dx**2)*lambdas_h[i+1]
                b = Cs-c-a

                if i == 0:
                    M[0][0] = a*(2*dx*alpha/lambdas_h[0])+b
                    M[0][1] = a + c

                elif i == Nx-1:
                    M[end_index][end_index-1] = a
                    M[end_index][end_index] = b + c
                else:
                    M[i][i-1] = a
                    M[i][i] = b
                    M[i][i+1] = c

        Initial = np.copy(sp.sparse.csr_matrix(sp.linalg.inv(M)).dot(D))


    for i in range(Nx):
        sheet.cell(row=i+1, column=1).value = float(Initial[i]) - 273
        sheet.cell(row=i+1, column=2).value = L - i*dx

    book.save("Numerical_solution.xlsx")
    book.close()

numerical_solution()

# def main():

