using Plots
using LinearAlgebra

Δx = 0.1
Δy = 0.1
x = Array(Δx:Δx:1-Δx) # first and last point already specified, 0)
y = Array(Δy:Δy:1-Δy) # first and last point already specified, 0
N = length(x)

f(x,y) = sin(2π*x) * y^2

function buildGrid(x,y)
    dimA = length(x), length(y)
    A = zeros(dimA[1], dimA[2])
    for ix in 1:1:dimA[1], iy in 1:1:dimA[2]
        A[ix, iy] = [ix*Δx, iy*Δy]
    end
    print(A)
    return A
end


buildGrid(x,y)
