module Chebyshev
using LinearAlgebra
export chebGrid, chebGridPoint, chebT, chebCardinal, coeffsToFunction, derivative

function chebGridPoint(i::Int, N::Int)
    return cos(Float64(i)*π / Float64(N))
end

"""grid of N+1 points"""
function chebGrid(N::Int)
    i = 0:N
    return chebGridPoint.(i, N)
end

function chebT(n::Int)
    if n < 0
        error("n must be >= 0")
    elseif n == 0
        return x -> 1
    elseif n == 1
        return x -> x
    else
        return x -> 2*x*chebT(n-1)(x) - chebT(n-2)(x)
    end
end

"""p(0)=p(N)=2, p(j)=1"""
function pj(j::Int, N::Int)
    j < 0 && error("j must be >=0")
    j > N && error("j must be <=N")
    j == 0 && return 2
    j == N && return 2
    return 1
end

function chebCardinalSummand(j::Int, m::Int, N::Int, grid::Array)
    ret = x -> 1.0/Float64(pj(m, N)) * chebT(m)(grid[j+1]) * chebT(m)(x)
    return ret
end

function chebCardinal(j::Int, grid::Array)
    N = length(grid) - 1
    factor = 2.0 / Float64(N * pj(j, N))
    foo = chebCardinalSummand(j, 0, N, grid)
    for m = 1:N
        foo = let fooOld = foo
            x -> fooOld(x) + chebCardinalSummand(j, m, N, grid)(x)
              end
    end
    return x -> factor * foo(x)
end

"""from Chebyshev coefficients return function that can
be called for any argument (in intervall -1..1)"""
function coeffsToFunction(coeffs::Array, grid::Array)
    N = length(grid) - 1
    foo = x -> coeffs[0 + 1] * chebCardinal(0, grid)(x)
    for j = 1:N
        foo = let fooOld = foo
            x -> fooOld(x) +  coeffs[j + 1] * chebCardinal(j, grid)(x)
        end
    end
    return foo
end

function dijValues(i::Int, j::Int, grid::Array)
    N = length(grid) - 1
    i == 1 && j == 1 && return (1 + Float64(2*(N^2)))/6.0
    i == N + 1 && j == N + 1 && return -(1 + Float64(2*(N^2)))/6.0
    i == j && j > 1 && j < N + 1 && return -grid[j] / (2.0 * (1.0 - grid[j]^2))
    return (-1)^(i + j - 2) * pj(i - 1, N) / (pj(j - 1, N) * (grid[i] - grid[j]))
end

function derivative(grid; k::Int = 1)
    N = length(grid) - 1
    Dij = Array{Float64,2}(undef, N+1, N+1)
    for ij in CartesianIndices(Dij)
        i = ij[1]; j = ij[2]
        Dij[i, j] = dijValues(i, j, grid)
    end
    return Dij^k
end

end
