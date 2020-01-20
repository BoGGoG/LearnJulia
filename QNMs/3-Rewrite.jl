using LinearAlgebra, Plots, ApproxFun
include("Chebyshev2.jl")

space = Chebyshev(0..1)

##
# equation:
# c00 has no w dependence
# c01 has the w dependence
# c00(u) ϕ(u) + c01(u,w) + c10(u) ϕ'(u) + c11(u,w) ϕ'(u) + c20(u) ϕ''(u) + c21(u,w) ϕ''(u)

Dij = Derivative(space)
DDij = Derivative(space, 2)
dspace = rangespace(Dij) # derivative space (ultraspherical)
ddspace = rangespace(DDij) # derivative space 2 (ultraspherical 2)
Sstd = Conversion(space, dspace) # space to derivative space
Sstdd = Conversion(space, ddspace) # space to second derivative space
Sdtdd = Conversion(dspace, ddspace) # derivative space to second

u = Fun(space)
q = 0.0

c00 = Multiplication(-3 - 9u^2 - 4q^2 * u^2, space)
c01 = Multiplication(6im*u, space)
c10 = Multiplication(u*(3 - 7u^4),
c11 = 4im*u^2
c20 = u*(u - u^5)
c21 = 0
