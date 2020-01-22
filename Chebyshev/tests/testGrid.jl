using Test

@testset "test hello world" begin
    @test helloWorld() == 1
end

@testset "test grid" begin
    N = 10 # number of grid points
    i = 3 # some grid point
    x = 0.3 # some point
    @test chebGridPoint(i, N) ≈ cos(i * π / N)
    grid = chebGrid(N)
    @test length(grid) == N + 1
end


# Notes
# =============
# - isapprox(a, b, atol = 1e-5)
