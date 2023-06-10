#!/usr/bin/env python3

mu = 1.983e-5
kappa = 0.02624
cp = 1006
g = 9.81
rho = 1.007
beta = 0.0037
area = 9.23*2.46
perimeter = 2 * 9.23 + 2*2.46
delta_t = 4

L = 2.46  # area/perimeter
Pr = mu * cp/kappa
Gr = L * L * L * rho * rho * g * delta_t * beta / (mu * mu)
Ra = Gr * Pr

t1 = (0.492/Pr)**(9.0/16.0)
t2 = (1 + t1)**(8.0/27)
t3 = 0.387 * Ra**(1.0/6.0)
t4 = 0.825 + t3 / t2

Nu = t4 * t4

h = Nu * kappa / L

print("Pr = ", Pr)
print("Gr = ", Gr)
print("Ra = ", Ra)
print("Nu = ", Nu)
print("h  = ", h)
