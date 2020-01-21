module ChebyshevTest
include("Chebyshev.jl")
using .Chebyshev
using Plots

allTests = false
# allTests = true

if allTests
    ## chebGridPoint
    chebGridPoint(12, 20)

    ## chebGrid
    g = chebGrid(30)
    print(typeof(g))
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
    # display(pl)

    ## pj
    println(Chebyshev.pj(0, 10) == 2)
    println(Chebyshev.pj(10, 10) == 2)
    println(Chebyshev.pj(1, 10) == 1)

    ## chebCardinal
    #
    N = 15
    grid = chebGrid(N)
    println(chebCardinal(1, grid)(0.3))

    xPoints = -1:0.08:1
    pl = plot()
    for i in 0:6
        plot!(pl, xPoints, chebCardinal(i, grid).(xPoints))
    end
    display(pl)

    ## Function as sum over cardinal functions
    fun(x) = 1 + x^2 * cos(x)
    pl = plot(xPoints, fun.(xPoints), label = "function")
    display(pl)

    funGrid = fun.(grid)
    scatter!(pl, grid, funGrid, label = "grid points")
    display(pl)

    scatter!(pl, xPoints, coeffsToFunction(funGrid, grid).(xPoints), label = "Cheb Interpolation")
    display(pl)
end



# N = 12
# grid = chebGrid(N)
# xPoints = -1:0.02:1
# fun(x) = cos(10*cos(x^2))

# pl = plot(xPoints, fun.(xPoints), label = "function")
# display(pl)

# coefs = fun.(grid)
# scatter!(pl, grid, coefs, label = "grid points")
# display(pl)

# scatter!(pl, xPoints, coeffsToFunction(coefs, grid).(xPoints), label = "Cheb Interpolation")
# display(pl)

## Derivatives

using ForwardDiff

N = 15
grid = chebGrid(N)
xPoints = -1:0.03:1
f(x) = cos(10*cos(x^2))
coeffs = f.(grid)
# plt = scatter(xPoints, f.(xPoints), label = "function")

fCheb = coeffsToFunction(coeffs, grid)
df(x) = ForwardDiff.derivative(f, x)

plt = plot(f, -1, 1, label = "f(x)")
plot!(plt, fCheb, label = "Chebyshev approximation")


D = derivative(grid)
dcoeffs = D*coeffs
dfCheb = coeffsToFunction(dcoeffs, grid)
plt2 = plot(df, label = "f'(x)", -1, 1)
plot!(plt2, dfCheb, label = "f'(x) Chebyshev approximation")

# plt = plot()
# # plt = scatter!(plt, plotPoints, fCheb.(plotPoints), label = "Cheb Interpolation")
# plot!(plt, f, label = "f(x)")
# plot!(plt, df, label = "f'(x)")
# plot!(plt, plotPoints, dfCheb.(plotPoints), label = "Cheb Interpolation")

display(plt)
display(plt2)



end
