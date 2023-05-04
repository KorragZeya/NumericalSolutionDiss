import numpy as np
import scipy as sp
import openpyxl

book = openpyxl.Workbook()

sheet = book.active

def numerical_solution():
    L = 3
    Nx = 100
    end_index = Nx - 1
    t = 20 * 86000
    dt = 100
    Nt = int(t/dt)
    dx = L/Nx
    Initial = np.full(Nx, 283)
    lambdas = np.full(Nx+1, 2.7)
    Cs = 2006200
    Cw = 4200000
    F = 0
    q = lambda i: 120*np.sin((2*3.14/86400)*i*dt)
    M = np.zeros((Nx, Nx), dtype=float)


    for k in range(Nt):
        D = np.array([Cs*Initial[0] + (dt/dx)*q(k) if i == 0
                      else Cs*Initial[i] for i in range(Nx)], dtype=float)

        for i in range(Nx):
                a = -(dt/dx**2)*lambdas[i]
                c = -(dt/dx**2)*lambdas[i+1]
                b = Cs-c-a
                if i == 0:
                    M[0][0] = a+b
                    M[0][1] = c
                elif i == Nx-1:
                    M[end_index][end_index-1] = a
                    M[end_index][end_index] = b + c
                else:
                    M[i][i-1] = a
                    M[i][i] = b
                    M[i][i+1] = c

        Initial = np.copy(sp.sparse.csr_matrix(sp.linalg.inv(M)).dot(D))

    print(type(float(Initial[2])))

    for i in range(Nx):
        sheet.cell(row=i+1, column=1).value = float(Initial[i])
        sheet.cell(row=i+1, column=2).value = L - i*dx

    book.save("Numerical_solution.xlsx")
    book.close()

numerical_solution()

# def main():

