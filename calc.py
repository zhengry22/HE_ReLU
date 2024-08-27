import numpy as np

# Define the polynomial coefficients
coeffs = [0.0020541, 0, -0.0207295, 0, 0.249792, 0.5, 0]

# Construct the polynomial
poly = np.poly1d(coeffs)

# Calculate f(x) for a specific value of x, for example, x = 2
x = float(input())
result = poly(x)

print(f"f({x}) = {result}")
