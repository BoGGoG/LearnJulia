using LinearAlgebra, Plots, ApproxFun
include("QNMs/Chebyshev2.jl")

space = Chebyshev(0..1)
f = Fun(x -> cos(10cos(x^2)), space)

## Check deriatives
plot(ylims = [-5, 10])
plot!(f)
plot!(f')
Dij = Derivative(space)
DDij = Derivative(space, 2)
plot!(Dij*f, linestyle = :dot, linewidth = 3 )
plot!(f'')
plot!(Dij*Dij*f, linestyle = :dot, linewidth = 3)
plot!(DDij*f, linestyle = :dash, linewidth = 2)
##


## Compare to analytic result (autodiff)
using ForwardDiff
foo(x) = cos(10*cos(x^2))
xpoints = 0:0.01:1
fooDiv(x) = ForwardDiff.derivative(foo, x)
fooDiv(0.5)
# plot(foo, xlims = [0,1])
plot(ylims = [-5,10])
# plot!(f')
plot!(Dij*f, linewidth =1, label = "D*f (Cheby)")
plot!(fooDiv, label = "autodiff", linestyle = :dot, linewidth = 4)


## Compare function with different orders of the polynomials (grid sizes)
f = Fun(x -> cos(10cos(x^2)), space)
f.coefficients
p = plot(f)
for i in 3:1:8
    g = Fun(space, f.coefficients[1:i])
    plot!(g)
end

## Derivative matrix acting on coefficients
# note that we have to use the Ultraspherical space as base for the
# function. I think it should consist of the derivatives of the
# Chebyshev polynomials
p2 = plot(f', label = "f'", legend = :bottomleft)
derivSpace = rangespace(Dij)
for i in 5:2:10
    h = Fun(derivSpace, Dij[1:i, 1:i] * f.coefficients[1:i])
    plot!(p2, h, label = "n = $i")
end
p2

## second derivative
p3 = plot(f'', label = "f''", legend = :bottomleft)
dderivSpace = rangespace(DDij)
for i in 5:5:10
    h = Fun(dderivSpace, DDij[1:i, 1:i] * f.coefficients[1:i])
    plot!(p3, h, label = "DDij, n = $i")
    g = Fun(dderivSpace, Dij[1:i, 1:i] * Dij[1:i, 1:i] * f.coefficients[1:i])
    plot!(p3, g, label = "Dij*Dij, n = $i")
end
p3

Dij

print(Derivative(space))
