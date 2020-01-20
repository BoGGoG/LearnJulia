module Chebyshev
using LinearAlgebra
export chebGrid, chebGridPoint, chebT, chebCardinal, coeffsToFunction

function chebGridPoint(i::Int, N::Int)
    return cos(Float64(i)*π / Float64(N))
end

"""grid of N+1 points"""
function chebGrid(N::Int)
    i = 0:N
    return chebGridPoint.(i, N)
end

function chebT(n::Int)
    # n < 0 && error("n must be >= 0")
    # n == 0 && return returnOne
    # n == 1 && return returnX
    # return 2*returnX*chebT(n-1) - chebT(n-2)
    #

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

end
