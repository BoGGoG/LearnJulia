module ChebyshevTest
include("Chebyshev.jl")
using .Chebyshev
using Plots

## chebGridPoint
chebGridPoint(12, 20)

## chebGrid
g = chebGrid(30)
pl = plot(g)
# display(pl)

## Chebyshev Polynomials
println(chebT(0)(123.0))
println(chebT(1)(2.0))
foo = chebT(2)
println(foo)
println(foo(2))
xPoints = -1:0.01:1
pl = plot()
for i in 1:6
    plot!(pl, xPoints, chebT(i).(xPoints))
end
display(pl)


end
