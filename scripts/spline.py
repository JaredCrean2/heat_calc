#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('AGG')

T_lower = -1
T_upper = 1

f_upper = 1
fprime_upper = 1

f_lower = -1
fprime_lower = 1

T_vals = np.linspace(-2, 2, 50);
f_vals = np.zeros(T_vals.shape)

a = f_lower
b = fprime_lower
c = 3*f_upper - 3*a - 2*b - fprime_upper  #3*(f_upper - f_lower) - 2*fprime_lower - fprime_upper
d = (fprime_upper - b - 2*c)/3  #2*(f_lower - f_upper) + fprime_lower + fprime_upper

a = f_lower
b = fprime_lower
c = 3*(f_upper - f_lower)  - 2*fprime_lower - fprime_upper;
d = fprime_upper + fprime_lower - 2*(f_upper - f_lower);

print ("from algebra, solution = ", [a, b, c, d])

#A = np.zeros((4, 4));
#b = np.zeros(4);
#
#A[0, 0] = 1;
#A[1, 1] = 1;
#A[2, :] = [1, 1, 1, 1];
#A[3, :] = [0, 1, 2, 3];
#
#b[0] = f_lower
#b[1] = fprime_lower
#b[2] = f_upper
#b[3] = fprime_upper;
#
#x = np.linalg.solve(A, b)
#print("from matrix, x = ", x)

for i in range(len(T_vals)):
  t = (T_vals[i] - T_lower)/(T_upper - T_lower)
  f_vals[i] = a + b*t + c*t*t + d*t*t*t
  print("i = ", i, ", t = ", t, ", f = ", f_vals[i])



fig, ax = plt.subplots(1);
fig.set_size_inches(8, 6)

ax.plot(T_vals, f_vals)
ax.grid()
fig.savefig("spline.png", dpi=300);



