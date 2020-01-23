using Test

@testset "test grid and functions" begin
    N = 15 # number of grid points
    i = 3 # some grid point
    x = 0.3 # some point
    grid = chebGrid(N)
    f(x) = cos(10*cos(x^2))
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
        coeffs = f.(grid)
        chebApproxFun = coeffsToFunction(coeffs, grid)
        @test f.(grid) ≈ chebApproxFun.(grid)
    end
end



# Notes
# =============
# - isapprox(a, b, atol = 1e-5)
# - test_trowhs() when error expected
