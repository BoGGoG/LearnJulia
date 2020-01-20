module Chebyshev
using LinearAlgebra
export chebGrid, chebGridPoint, chebT

function chebGridPoint(i::Int, N::Int)
    return cos(Float64(i)*π / Float64(N))
end

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




end
