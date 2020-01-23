using Test
using ForwardDiff

@testset "test derivatives" begin
    N = 20 # number of grid points
    i = 3 # some grid point
    x = 0.4 # some point
    f(x) = cos(10*cos(x^2))
    grid = chebGrid(N)
    coeffs = f.(grid)
    D = derivative(grid)

    # check if the derivative matrix has correct size,
    # take derivative of a function and check with
    # automatic differentiation of that function
    @testset "derivative matrix" begin
        @test size(D) == (N+1, N+1)
        d = size(coeffs)[1]
        @test size(D) == (d,d)
    end

    @testset "derivative function" begin
        fCheb = coeffsToFunction(coeffs, grid)
        df(x) = ForwardDiff.derivative(f, x)

        dcoeffs = D*coeffs
        dfCheb = coeffsToFunction(dcoeffs, grid)

        rtol = 0.01 # actually pretty bad
        @test isapprox(df(x), dfCheb(x), rtol = rtol)
    end
end
