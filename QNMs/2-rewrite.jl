using LinearAlgebra, Plots, ApproxFun
using StatsBase
include("Chebyshev2.jl")

function getQNMs(nGrid; space = Chebyshev(Domain(0..1)), q = 0.0)
    Dij = Derivative(space)
    DDij = Derivative(space, 2)
    u = Fun(space)
    q = 0.0

    dspace = rangespace(Dij) # derivative space (ultraspherical)
    ddspace = rangespace(DDij) # derivative space 2 (ultraspherical 2)
    Sstd = Conversion(space, dspace) # space to derivative space
    Sstdd = Conversion(space, ddspace) # space to second derivative space
    Sdtdd = Conversion(dspace, ddspace) # derivative space to second

    c00 = Multiplication(-3 - 9u^2 - 4q^2 * u^2, space)
    c01 = Multiplication(6im*u, space)
    c10 = Multiplication(u*(3 - 7u^4), dspace)
    c11 = Multiplication(4im*u^2, dspace)
    c20 = Multiplication(u*(u - u^5), ddspace)
    c21 = Multiplication(0*u, ddspace)

    # convert all to Ultraspherical(2, domain)
    M0 = Sstdd*c00 + Sdtdd*c10*Dij + c20*DDij
    M1 = Sstdd*c01 + Sdtdd*c11*Dij + c21*DDij
    λ = eigvals(Matrix(M1[1:nGrid, 1:nGrid]), -Matrix(M0[1:nGrid, 1:nGrid]))
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

qnms = getQNMs.(40:5:45)
qnms = removeSpuriousQNMs(qnms; digits = 4)
scatter(2π*qnms)
print(2π*qnms)
print("------------------")


qnms = getQNMs(300)
scatter(qnms)
