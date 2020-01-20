using LinearAlgebra, SpecialFunctions, Plots, ApproxFun

x = Fun(identity, 0..10)
f = sin(x^2)
g = cos(x)

h = f + g^2
r = roots(h)
rp = roots(h')


plot(h, label = "f + g^2")
scatter!(r, h.(r); label = "roots")
scatter!(rp, h.(rp); label = "extrema")

plot(h - h')
norm(h - h') # what is the norm of a function?

(h - h')(10)

g = cumsum(f)
D = Derivative()
B = Dirichlet

f2 = D*h
f3 = h'
plot(f2 - f3)

plot(f2)
plot!(f3)

x = Fun(identity, -1000..200)
D = Derivative()
B = Dirichlet()
L = D^2 - x

u = [B;L] \ [[1,2], 0]
plot(u; label = "u")

a = Fun(cos, Chebyshev())
b = Fun(x -> cos(10cos(x^2)), Chebyshev())
ncoefficients(b)
ncoefficients(a*b)

d = Domain(-1..1)^2
x,y = Fun(identity, d)
x = Fun(identity, 0..10)
f = exp.(-10(x+0.3)^2-20(y-0.2)^2)  # use broadcasting as exp(f) not implemented in 2D
sin(x)
Δ = Laplacian(d)
A = [Dirichlet(d);Δ]              # Δ is an alias for Laplacian()
@time u = A \ [zeros(∂(d));f]     #4s for ~3k coefficients

println(d)

c = 0.
x = Fun()
u0 = 0*x
ε = 0.1
N = u -> [u(-1.)-c; u(1.); ε*u'' + 6*(1-x^2)*u' + u^2 - 1.0]
u = newton(N, u0)
plot(u)
