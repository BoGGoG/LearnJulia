# following [1] to calculate solutions of first order differential equations
# using their Chebyshev method
# [1] https://arxiv.org/pdf/1202.1347.pdf

## 1) u'(x) + x^3 u(x) = 100sin(20000x^2)
# u(-1) = 0
space = Chebyshev(-1..1)
x = Fun(identity, space)
D = Derivative(space)
B = Dirichlet(space, 0)
a = Fun(t -> 100sin(20000t^2), space)
L = D + x^3
uCheb = [B;L] \ [0.; a]
plot(uCheb)


## 2) u''(x) + x^2/2 u(x) = 0
# u(-8) = u(8) = 0
s = Chebyshev(-8..8)
x = Fun(identity, s)
L = -Derivative(s)^2 + x^2/2
B = Dirichlet(s)
λ, v = ApproxFun.eigs(B, L, 500, tolerance = 1E-10)
λ

plot(v[1:3])

vCheb = [B;L] \ [0.; 0]
plot(vCheb)

## 3) 2, but on my own
