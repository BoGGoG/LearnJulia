using Plots
using LinearAlgebra

Δx = 0.1
x = Δx:Δx:1-Δx # first and last point already specified, 0
N = length(x)
B = sin.(2π * x)
A1 = zeros(N, N)

function buildLaplace1(N::Int)
    A = zeros(N, N)
    for i = 1:N, j = 1:N
        abs(i - j) <= 1 && (A[i, j] += 1)
        i == j && (A[i, j] -= 3)
    end
    return A
end

function buildLaplace2(N::Int)
    A = zeros(N, N)
    for i = 1:N, j = 1:N
        abs(i - j) == 1 ? A[i, j] = 1 : 1
        i == j ? A[i, j] = 2 : 1
    end
    return A
end

function buildLaplace3(N::Int)
    diagonal = -2 * ones(N)
    diagonal2 = 1 * ones(N - 1)
    A = Tridiagonal(diagonal2, diagonal, diagonal2)
    return A
end


function buildLaplace4(N::Int)
    A = SymTridiagonal(fill(-2, N), fill(1, N - 1))
    return A
end

println("Compiling")
@time buildLaplace1(5)
@time buildLaplace2(5)
@time buildLaplace3(5)
@time buildLaplace4(5)

println("---------------------")
println("Testing with N=10000")
@time buildLaplace1(10000)
@time buildLaplace2(10000)
@time buildLaplace3(10000)
@time buildLaplace4(10000)

A = buildLaplace4(N)
A = A/(Δx^2)

U = A\B
plot([0;x;1], [0;U;0], label="U")
plot!([0;x;1], -sin.(2π*[0;x;1])/(4π^2), label = "Analytical Solution")


# finite difference method
using DiffEqOperators
using BandedMatrices
order = 2
deriv = 2
Δx = 0.1
N = 9
A = CenteredDifference(deriv, order, Δx, N)
printstyled(BandedMatrix(CenteredDifference(4,2,Δx,N)))

using Pkg
pkg"add FFTW"
using FFTW
x = range(0, 2π, length = 100)
u(x) = sin(x)
freqs = fft(u.(x))[1:length(x)÷2 + 1]
c = 2*abs.(freqs/length(x))'
freqRange = 0:(length(c) - 1)
scatter(freqRange, c')

function buildFFT(c, N, x)
    fx = zeros(length(x))
    for i in range(0, stop = (N-1))
        fx .+= c[i+1]*sin.(i*x)
    end
    return fx
end

y = buildFFT(c, 4, x)
scatter(x,y)

using ApproxFun
S = Fourier()
n = 100
T = ApproxFun.plan_transform(S, n)
Ti = ApproxFun.plan_itransform(S, n)
x = points(S, n)

r = (T*cos.(cos.(x.-0.1)))
# r = (cos.(cos.(x.-0.1)))
scatter(x,r)
y = buildFFT(r, 5, x)
scatter(x,y)
scatter(x, r[1].*sin.(0*x) + r[2].*sin.(1*x) + r[3].*sin.(2*x) + r[4].*sin.(3*x))

D2 = Derivative(S, 2)
D2

# Chebyshev
S = Chebyshev()
n = 100
T = ApproxFun.plan_transform(S, n)
Ti = ApproxFun.plan_itransform(S, n)
x = points(S, n)
r = (T*cos.(cos.(x.-0.1)))
scatter(x,r)
D2 = Derivative(S,2);

println(r)
println(x)


# AU=B problem
using Pkg
pkg"add IterativeSolvers"
n = 10
A = rand(n,n)
B = rand(n)
using IterativeSolvers, LinearAlgebra
U = gmres(A, B, tol = 1e-8)
norm(A*U - B)
A*U - B

# Poinsson Equation
struct SizedStrangMatrix
    size::Tuple{Int,Int}
end
Base.eltype(A::SizedStrangMatrix) = Float64
Base.size(A::SizedStrangMatrix) = A.size
Base.size(A::SizedStrangMatrix, i::Int) = A.size[i]

A = SizedStrangMatrix((length(B), length(B)))

import LinearAlgebra
function LinearAlgebra.mul!(C, A::SizedStrangMatrix, B)
    for i in 2:length(B)-1
        C[i] = B[i-1] - 2B[i] + B[i+1]
    end
    C[1] = -2B[1] + B[2]
    C[end] = B[end-1] - 2B[end]
    C
end

Base.:*(A::SizedStrangMatrix, B::AbstractVector) = (C = similar(B); mul!(C,A,B))

using IterativeSolvers, LinearAlgebra
U = gmres(A,B, tol  = 1e-14)
norm(A*U - B)

# Spectal time stepping for heat equation in fourier space
using Pkg
pkg"add ApproxFun Sundials Plots DifferentialEquations"

using ApproxFun, Sundials, Plots; gr()

S = Fourier()
n = 100
x = points(S, n)
T = ApproxFun.plan_transform(S, n)
Ti = ApproxFun.plan_itransform(S, n)

# convert inital condition to Fourier Space
u0 = T*cos.(cos.(x.-0.1))
D2 = Derivative(S, 2)
L = D2[1:n, 1:n]

using LinearAlgebra
heat(du, u, L, t) = mul!(du, L, u)

prob = ODEProblem(heat, u0, (0.0, 10.0), L)

sol = solve(prob, CVODE_BDF(linear_solver=:Diagonal); reltol = 1e-8, abstol = 1e-8)

plot(x, Ti*sol(0.0))
plot!(x, Ti*sol(0.5))
plot!(x, Ti*sol(2.0))
plot!(x, Ti*sol(10.0))
