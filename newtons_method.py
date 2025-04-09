import numpy as np

def f(x):
    # Define the function for which we want to find the root
    return np.log(np.exp(x) + np.exp(-x))

def f_prime(x):
    # Define the derivative of the function
    return (np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x))

def newtons_method(f, f_prime, x0, tolerance=1e-7, max_iterations=1000):
    x = x0
    for _ in range(max_iterations):
        x_new = x - f(x) / f_prime(x)
        if abs(x_new - x) < tolerance:
            return x_new
        x = x_new
        print(x)
    raise ValueError("Newton's method did not converge")

# Example usage
initial_guess = 1.0
root = newtons_method(f, f_prime, initial_guess)
print(f"The root is: {root}")