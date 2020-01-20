using LinearAlgebra, SpecialFunctions, Plots, ApproxFun
using GenericLinearAlgebra
using StatsBase

# space = Chebyshev()
space = Chebyshev(Domain(0..1))
Dij = Derivative(space)
DDij = Derivative(space, 2)
u = Fun(space)
q = 0.0

c00 = -3 - 9u^2 - 4q^2 * u^2
c01 = 6im*u
c10 = u*(3 - 7u^4)
c11 = 4im*u^2
c20 = u*(u - u^5)
c21 = 0

M0 = c00*I + c10*Dij + c20*DDij
M1 = c01*I + c11*Dij + c21*DDij
# M = M0 + ω*M1

# p = plot(xlims = (-15,15), ylims = (-10,10))
p = plot(alpha = 0.1)
for n in 200:20:300
    # λ, v = eigen(Matrix(M1[1:n, 1:n]), Matrix(M0[1:n, 1:n]))
    λ = eigvals(Matrix(M1[1:n, 1:n]), Matrix(M0[1:n, 1:n]))
    println(4π*λ[end-5:end])
    println("-------")
    scatter!(p, 2π*λ)
end
p
λ[end-3:end]

num = 100
λ, v = eigen(Matrix(M0[1:num, 1:num]), Matrix(M1[1:num, 1:num]))


big.(M1[1:100, 1:100])

eigvals(big.(randn(4,4)))

function getQNMs(nGrid; space = Chebyshev(Domain(0..1)), q = 0.0)
    Dij = Derivative(space)
    DDij = Derivative(space, 2)
    u = Fun(space)

    c00 = -3 - 9u^2 - 4q^2 * u^2
    c01 = 6im*u
    c10 = u*(3 - 7u^4)
    c11 = 4im*u^2
    c20 = u*(u - u^5)
    c21 = 0

    M0 = c00*I + c10*Dij + c20*DDij
    M1 = c01*I + c11*Dij + c21*DDij
    λ = eigvals(Matrix(M1[1:nGrid, 1:nGrid]), Matrix(M0[1:nGrid, 1:nGrid]))
    return λ
end

function removeSpuriousQNMs(qnms::Vector{Array{Complex{Float64}, 1}}; digits = 3)
    qnms = copy(qnms)
    qnms = map(vec -> round.(vec, digits = digits), qnms)
    qnms = collect(Iterators.flatten(qnms))
    qnmsSorted = sort(qnms, lt = (x,y) -> real(x) < real(y))
    duplicates = []
    nonunique(x) = [k for (k, v) in countmap(x) if v > 1]
    return nonunique(qnmsSorted)
end

qnms = getQNMs.(800:10:810)
qnms = removeSpuriousQNMs(qnms; digits = 6)
scatter(4π*qnms)
