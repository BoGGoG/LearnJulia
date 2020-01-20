using Plots
using LinearAlgebra

Δx = 0.1
x = Δx:Δx:1-Δx # first and last point already specified, 0
N = length(x)
B = sin.(2π * x)
A1 = zeros(N, N)

function buildLaplace(N::Int)
    A = SymTridiagonal(fill(-2, N), fill(1, N - 1))
    return A / (Δx^2)
end

Δ = buildLaplace(N)

# Δu = b
# b = Δ\b
u = Δ\B

plot(u)
plot!(-B ./ (4 * π^2))
