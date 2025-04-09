import numpy as np

def dispersion_free_surface(alpha, N, h=1):
    """
    Calculates the positive imaginary and N first positive real solutions of -alpha = k*tan(k h)
    for complex alpha. It uses three methods - homotopy for starting with alpha =1, guess with linear expansion
    and a linear expansion. The first roots is positive imaginary and the next are the first N positive real
    ordered from smallest.
    
    If the value for h is not given the default value is h = 1.
    """
    if h != 1:
        alpha = h * alpha  # scaling for the solution

    mroots = np.zeros(N + 1, dtype=complex)

    if N == 0:  # the N = 0 case does not involve any of the special methods and is treated separately.
        count = 0
        mroots[count] = homotopy(alpha, count)
    else:
        count = 0
        mroots[count] = homotopy(alpha, count)
        count += 1
        while True:
            mroots[count] = homotopy(alpha, count)
            if abs(mroots[count] - (1j * count * np.pi + alpha / (1j * count * np.pi))) < 0.01:
                while True:
                    mroots[count] = oneroot(alpha, 1j * count * np.pi + alpha / (1j * count * np.pi))
                    if abs(mroots[count] - (1j * count * np.pi + alpha / (1j * count * np.pi))) < 1e-8:
                        mroots[count:N + 1] = 1j * np.arange(count, N + 1) * np.pi + alpha / (1j * np.arange(count, N + 1) * np.pi)
                        count = N
                        break
                    if count == N:
                        break
                    count += 1
            if count == N:
                break
            count += 1

    mroots = 1 / h * mroots
    mroots[0] = mroots[0]
    return mroots

def homotopy(alpha, N):
    """
    Calculates the Nth root using the homotopy method.
    """
    if N == 0:
        mroot = oneroot(1, 1)
    else:
        mroot = oneroot(1, 1j * N * np.pi)

    step = 0.043
    if abs(alpha) < 1:
        alphastep = np.concatenate(([1], np.arange(1 - step, abs(alpha), -step), [abs(alpha)]))
    else:
        alphastep = np.concatenate(([1], np.arange(1 + step, abs(alpha), step), [abs(alpha)]))

    for k in range(1, len(alphastep)):
        mroot = oneroot(alphastep[k], mroot)

    if np.angle(alpha) > 0:
        alphastep = abs(alpha) * np.exp(1j * np.arange(0, np.angle(alpha) + np.pi / 30, np.pi / 30))
    else:
        alphastep = abs(alpha) * np.exp(1j * np.arange(0, np.angle(alpha) - np.pi / 30, -np.pi / 30))

    for k in range(1, len(alphastep)):
        mroot = oneroot(alphastep[k], mroot)

    return mroot

def oneroot(alpha, guess):
    """
    Calculates the root nearest the root guess.
    """
    ans1 = guess + 1
    out = guess
    while abs(ans1 - out) > 1e-9:
        ans1 = out
        out = ans1 - f(ans1, alpha) / difff(ans1)
    return out

def f(z, alpha):
    return z * np.tanh(z) - alpha

def difff(z):
    return np.tanh(z) + z * (1 / np.cosh(z)) ** 2
