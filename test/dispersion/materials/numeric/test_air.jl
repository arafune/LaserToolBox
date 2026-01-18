using Test
using LaserToolBox

@testset "Air (Numeric)" begin
    λ = 0.8
    @test n.air(λ) ≈ 1.0002750477973053 atol = 1e-10

    @testset "Derivative Consistency" begin
        dn = n.air(λ; derivative = 1)
        @test dn < 0

        d2n = n.air(λ; derivative = 2)
        @test d2n > 0

        @test_throws ArgumentError("derivative must be ≥ 0") n.air(λ; derivative = -1)
    end

    λ = 0.01
    @test_throws ArgumentError n.air(λ)

    @testset "Vectorization" begin
        λs = [0.4, 0.8, 1.2]
        @test n.air.(λs) isa Vector{Float64}
        @test issorted(n.air.(λs), rev = true)
    end
end

