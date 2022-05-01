from sympy.solvers import solve
from sympy import *

r1, r2 = symbols('r1 r2')
L = Symbol('L')
g1 = 1 - L/r2
g2 = 1 - L/r1
g1g2 = g1 * g2

print(factor(g1g2))
print(expand(g1g2))

eqn = g1g2 - 1

print(solve(eqn, L))
print(solve(eqn.subs(r1, 1.4).subs(r2, 1.4), L))
print(solve(eqn.subs(r1, 1.0).subs(r2, 1.4), L))
print(solve(eqn.subs(r1, oo).subs(r2, 1.4), L))
