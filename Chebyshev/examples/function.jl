include("../Chebyshev.jl")
using .Chebyshev
using Plots
using ForwardDiff

nGrid = 15
grid = chebGrid(nGrid)

f(x) = cos(10*cos(x^2))
coeffs = f.(grid)
fCheb = coeffsToFunction(coeffs, grid)

df(x) = ForwardDiff.derivative(f, x)
D = Chebyshev.derivative(grid)
dcoeffs = D*coeffs
dfCheb = coeffsToFunction(dcoeffs, grid)

plt = plot(f, -1, 1, label = "f(x) = cos(10*cos(x^2))")
plot!(plt, fCheb, -1, 1, linestyle = :dash, linewidth = 4, label = "Chebyshev approx, N = $nGrid")
plot!(plt, x -> df(x)/4.0, -1, 1, label = "f'(x)/4.0 automatic differentiation")
plot!(plt, x -> dfCheb(x)/4.0, -1, 1, linestyle = :dash, linewidth = 4, label = "f'(x)/4.0 Chebyshev, N = $nGrid")
display(plt)
