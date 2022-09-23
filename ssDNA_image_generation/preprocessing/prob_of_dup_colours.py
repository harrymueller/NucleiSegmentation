import numpy as np
import math

# numpy and math shortcuts
ull = np.ulonglong
ld  = np.longdouble
ln  = np.log
e   = math.e
pi  = math.pi

# n & r
n   = ld(math.pow(128, 3))
r   = ld(3500)

# stirlings and gosper's approximations
def stirlings_approx_ln(s):
    return s * ln(s) - s

def gosper_formula(s):
    return 0.5 * ln((2*s + 1/3)*pi) + s*ln(s) - s

ps = stirlings_approx_ln(n) - stirlings_approx_ln(n - r) - r * ln(n)
pg = gosper_formula(n) - gosper_formula(n - r) - r * ln(n)

print("Using Stirling's Approximation")
print(1 - np.exp(ps))
print("Using Gosper's Formula")
print(1 - np.exp(pg))