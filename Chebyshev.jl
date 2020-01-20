using ApproxFun

S = Chebyshev()
n = 100
T = ApproxFun.plan_transform(S, n)
Ti = ApproxFun.plan_itransform(S, n)
x = points(S, n)

r = (T*cos.(cos.(x.-0.1)))
r2 = Fun(x->cos(cos(x-0.1)), S)
D2 = Derivative(S, 2)


plot(x, D2*r2)

plot(x, r2)
plot(x, r)


f(x) = cos(cos(x-0.1))
xlin = range(0,1; length = 100)
plot(xlin, f.(xlin))

# EXAMPLE: Reaction Heat Equation in Fourier Space
S = Fourier()
n = 100
x = points(S, n)
T = ApproxFun.plan_transform(S, n)
Ti = ApproxFun.plan_itransform(S, n)

# Convert the initial condition to Fourier space
u₀ = T*cos.(cos.(x.-0.1))
D2 = Derivative(S,2)
L = D2[1:n,1:n]

using DiffEqOperators, LinearAlgebra, DifferentialEquations

A = DiffEqArrayOperator(Diagonal(L))

function f(dû,û,tmp,t)
  # Transform u back to point-space
  mul!(tmp,Ti,û)
  # apply nonlinear function 0.75sqrt(u)-u in point-space
  @. tmp = 0.75sqrt(tmp) - tmp
  mul!(dû,T,tmp) # Transform back to Fourier space
end

# Define u' = Au + f
prob = SplitODEProblem(A, f, u₀, (0.0,10.0), similar(u₀));
prob = SplitODEProblem(A, f, u0, (0.0, 10.0), similar(u0))
sol = solve(prob, KenCarp4())
