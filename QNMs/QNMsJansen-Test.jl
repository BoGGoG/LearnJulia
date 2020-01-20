include("QNMsJansen.jl")
using .QNMs
using LinearAlgebra
using ApproxFun
using Plots
const N = 20

##
x = QNMs.Chebyshev(0..1)
ps = points(Chebyshev(), N)
xarr = QNMs.chebGrid(N; space = "banana")
plot(ps)
plot!(xarr)
# why not the same?
##


foo = QNMs.chebT.([1,2,3])

plot(QNMs.chebT(1).(-1:0.05:1))
plot!(QNMs.chebT(2).(-1:0.05:1))
plot!(QNMs.chebT(3).(-1:0.05:1))
plot!(QNMs.chebT(4).(-1:0.05:1))
plot!(QNMs.chebT(5).(-1:0.05:1))
plot!(QNMs.chebT(6).(-1:0.05:1))

##

QNMs.pj(-1, N)
QNMs.pj(0, N)
QNMs.pj(1, N)
QNMs.pj(N, N)
QNMs.pj(3, N)


##

foo = x -> chebT(0)(x) + chebT(1)(x)
plot(foo.(ps))

foo2 = x -> sum([chebT(1)(x), chebT(2)(x)])
foo2(2)
foo3 = x -> sum(chebT.([1,2,3])(x))
cheby = chebGrid(N; space = "Chebyshev")

chebArrFun(m) = x -> chebT(m)(cheby[m+1]) * chebT(m)(x) / pj(m, N)
asdf = chebArrFun.([i for i in 1:5])
asdf[1](1)

function sumLambda(fooArray, x)
        tmp = 0
        for i in eachindex(fooArray)
                tmp += fooArray[i](x)
        end
        return tmp
end


foo(1)
