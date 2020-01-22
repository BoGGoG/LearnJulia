using .Chebyshev
using Test

# include("testGrid.jl")

@testset "test grid and functions" begin
    N = 10 # number of grid points
    i = 3 # some grid point
    x = 0.3 # some point
    grid = chebGrid(N)
    @testset "grid" begin
        @test chebGridPoint(i, N) ≈ cos(i * π / N)
        @test length(grid) == N + 1
    end

    @testset "chebT" begin
        @test_throws ErrorException("j must be >=0") Chebyshev.pj(-1, N)
        @test_throws ErrorException("j must be <=N") Chebyshev.pj(N + 1, N)
        @test Chebyshev.pj(0, N) == 2
        @test Chebyshev.pj(N, N) == 2
        @test typeof(chebT(1)(0.3)) == Float64
    end

    # check if cardinal functions at least work
    # and that the chebyshev transformed function
    # gives the correct values at the grid points
    @testset "chebCardinal" begin
        @test typeof(chebCardinal(i, grid)(0.3)) == Float64
        fun(x) = cos(10*cos(x^2))
        coeffs = fun.(grid)
        chebApproxFun = coeffsToFunction(coeffs, grid)
        @test fun.(grid) ≈ chebApproxFun.(grid)
    end

    # check if the derivative matrix has correct size,
    # take derivative of a function and check with
    # automatic differentiation of that function
    @testset "derivatives" begin
    end
end

# @testset "current @test_throwsdev" begin
# end

