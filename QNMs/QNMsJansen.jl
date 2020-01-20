module QNMs
using ApproxFun
using Plots
using LinearAlgebra


export chebGrid, chebT, pj

function chebGrid(nPoints::Int; domain = -1..1, space = "Chebyshev")
    if space != "Chebyshev"
        print(space, " space not implemented, using Chebyshev grid")
        return chebGrid(nPoints)
    end
    i = 0:1:(nPoints - 1)
    xarr = cos.(i*π/(nPoints - 1))
    return xarr
end

function chebT(n)
    if n == 0
        return one
    elseif n == 1
        return identity
    elseif n < 0
        print("T_i does not exist for i<0, returning 0")
        return zero
    else
        return x-> 2*identity(x)*chebT(n-1)(x) - chebT(n-2)(x)
    end
end

function pj(i::Int, nPoints::Int)
    if i < 0
        print("Error in pj(i), does not exist for i<0. Returning 0")
        return 0
    elseif i==0 || i==nPoints
        return 2
    else
        return 1
    end
end


function sumLambda(fooArray, x)
        tmp = 0
        for i in eachindex(fooArray)
                tmp += fooArray[i](x)
        end
        return tmp
end

function cardinalFunction(j::Int, nPoints::Int)
    cheby = chebGrid(nPoints; space = "Chebyshev")
    f0 = x -> 2 / (float64(N)*pj(j,nPoints))
    f1 = x -> 0
    # for m in 0:N
    #     f1 = x -> f1(x)*1
    # end
    chebArrFun(m) = x -> 2 / (float64(N)*pj(j,nPoints)) * chebT(m)(cheby[m+1]) * chebT(m)(x) / pj(m, nPoints)
    chebArr = chebArrFun.([i for i in 1:j])
    f = x -> sum[chebArr(x)]
    return f
end



end
