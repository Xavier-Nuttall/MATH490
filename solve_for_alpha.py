import random
import numpy as np
from wave_numbers import dispersion_free_surface
from functools import lru_cache

import numpy as np

# Parameters
frequency = 0.001  # Example value for omega
gravity = 9.81  # Example value for g
alpha = frequency**2 / gravity
BARRIER_DEPTH = 10.0
MAX_DEPTH = 20.0
MAX_TRUNCATION = 250
STEP_TRUNCATION = 50
MAX_COLLOCATION_POINTS = 250
STEP_COLLOCATION_POINTS = 50

# Functions

@lru_cache(maxsize=None)
def phi_n_norm(n):
    return (2*k_n[n] * MAX_DEPTH + np.sinh(2*k_n[n]*MAX_DEPTH))/(4*k_n[n])


@lru_cache(maxsize=None)
def K_denominator(n):
    return phi_n_norm(n) * k_n[n] * 1j

@lru_cache(maxsize=None)
def phi_n(z, n):
    return np.cosh(k_n[n] * (z + MAX_DEPTH))

@lru_cache(maxsize=None)
def phi_0(z):
    return np.cosh(k_n[0] * (z + MAX_DEPTH))


def phiMinus(z, x, truncation=MAX_TRUNCATION):
    summation = sum(A_n[n] * np.exp(1j * k_n[n] * x) + B_n[n] * np.exp(-1j * k_n[n] * x) for n in range(truncation))
    return np.real(summation * np.cosh(k_n[0] * (z + MAX_DEPTH)))

def phiPlus(z, x, truncation=MAX_TRUNCATION):
    summation = sum(C_n[n] * np.exp(1j * k_n[n] * x) + D_n[n] * np.exp(-1j * k_n[n] * x) for n in range(truncation))
    return np.real(summation * np.cosh(k_n[0] * (z + MAX_DEPTH)))

# Solving for coefficients
def calculateB_n(n, integral_point_count):
    A = A_n[n]
    innerProduct = 0
    for i in range(integral_point_count):
        innerProduct += u_n[i] * phi_n(xi_points[i], n) * weights[i]

    return A - innerProduct/(phi_n_norm(n) * k_n[n] * 1j)

def calculateC_n(n, integral_point_count):
    D = D_n[n]
    innerProduct = 0
    for i in range(integral_point_count):
        innerProduct += u_n[i] * phi_n(xi_points[i], n) * weights[i]

    return D + innerProduct/(phi_n_norm(n) * k_n[n] * 1j)

# Integral kernel
K_cache = {}

def K(z, xi, truncation):
    if (z, xi) in K_cache:
        return K_cache[(z, xi)]
    if (xi, z) in K_cache:
        return K_cache[(xi, z)]
    
    output = np.complex128(0)
    for n in range(truncation):
        output += phi_n(z, n) * phi_n(xi, n) / K_denominator(n)
        
    K_cache[(z, xi)] = output
    return output

A_0 = 1

# Coefficients for waves
A_n = np.zeros(MAX_TRUNCATION, dtype=np.complex128)
D_n = np.zeros(MAX_TRUNCATION, dtype=np.complex128)
B_n = np.zeros(MAX_TRUNCATION, dtype=np.complex128)
C_n = np.zeros(MAX_TRUNCATION, dtype=np.complex128)


# Difference from actual solution at x = 0
def functionalDifference(truncation):
    return sum(A_n[:truncation]) + sum(B_n[:truncation]) - sum(C_n[:truncation]) - sum(D_n[:truncation])

A_n[0] = A_0



k_n = dispersion_free_surface(alpha, MAX_TRUNCATION, MAX_DEPTH) 



truncation_values = range(50, MAX_TRUNCATION, STEP_TRUNCATION)

collocation_points_range = range(50, MAX_COLLOCATION_POINTS, STEP_COLLOCATION_POINTS)

difference = np.zeros((len(truncation_values), len(collocation_points_range)), dtype=np.complex128)

for k,truncation in enumerate(truncation_values):
    for j,collocation_points in enumerate(collocation_points_range):
        print(f"Truncation Value: {truncation}, Colocation Points: {collocation_points}")
        # Define the interval and number of points
        integral_point_count = collocation_points

        z_points = np.linspace(-BARRIER_DEPTH, 0, collocation_points)
        xi_points = np.linspace(-BARRIER_DEPTH, 0, integral_point_count)
        weights = np.ones(integral_point_count) * (BARRIER_DEPTH / (integral_point_count - 1))
        weights[0] /= 2
        weights[-1] /= 2


        # Construct the matrix A and vector b
        A = np.zeros((collocation_points, integral_point_count), dtype=np.complex128)
        b = np.zeros(collocation_points, dtype=np.complex128)

        for i in range(collocation_points):
            b[i] = A_n[0] * phi_0(z_points[i])
            if i % 25 == 0:
                print(i)
            for n in range(integral_point_count):
                A[i, n] = K(z_points[i], xi_points[n], truncation)* weights[n]
                if np.imag(A[i, n]) > 1e-10:
                    print(f"Imaginary part of A[{i}, {n}] is {np.imag(A[i, n])}")

        


        # Solve the linear system
        print(A.shape, b.shape)
        u_n = np.linalg.solve(A, b)

        # Calculate B_n and C_n
        for n in range(truncation):
            B_n[n] = calculateB_n(n, integral_point_count)
            C_n[n] = calculateC_n(n, integral_point_count)

        difference[k, j] = functionalDifference(truncation)

        if k == 3 and j == 3:
            # print("A_n:", A_n)
            print("B_n:", B_n)
            # print("C_n:", C_n)
            # print("D_n:", D_n)
            # print("u_n:", u_n)
            print("Difference:", difference[k, j])
            
            print("Conservation of Energy:", B_n[0]**2 + C_n[0]**2, " == ", A_n[0]**2 + D_n[0]**2)
            # break
        

# Write difference to file
# np.savetxt("difference.txt", difference, delimiter=",")

# Plotting
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# print(difference)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# X, Y = np.meshgrid(truncation_values, collocation_points_range)
# ax.plot_surface(X, Y, np.abs(difference), cmap='viridis')
# ax.set_xlabel('Truncation Value')
# ax.set_ylabel('Collocation Points')
# ax.set_zlabel('Difference')
# plt.title('Difference in Solution')
# plt.show()
