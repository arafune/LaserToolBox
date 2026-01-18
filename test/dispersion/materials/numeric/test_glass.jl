using Test
using LaserToolBox

material_test_data = [
    (func = n.bk7, expected = 1.510776231419874, atol = 1e-10, λ = 0.8, label = "BK7"),
    (
        func = n.fused_silica,
        expected = 1.453317257044587,
        atol = 1e-10,
        λ = 0.8,
        label = "UVFS",
    ),
    (func = n.caf2, expected = 1.4305724647561817, atol = 1e-10, λ = 0.8, label = "CaF2"),
    (func = n.caf2, expected = 1.428, atol = 1e-3, λ = 1.064, label = "CaF2_1.064μm"),
    (func = n.sf10, expected = 1.7112, atol = 1e-4, λ = 0.800, label = "SF10"),
    (func = n.sf11, expected = 1.7646, atol = 1e-4, λ = 0.800, label = "SF11"),
    (func = n.sf11, expected = 1.7847, atol = 1e-4, λ = 0.5876, label = "SF11_dline"),
]

@testset "Glass Refractive Indices (Numeric)" begin
    for data in material_test_data
        @testset "$(data.label) (Numeric)" begin
            λ = data.λ
            @test data.func(λ) ≈ data.expected atol = data.atol

            @testset "Derivative Consistency" begin
                dn = data.func(λ; derivative = 1)
                @test dn < 0

                d2n = data.func(λ; derivative = 2)
                @test d2n > 0

                @test_throws ArgumentError("derivative must be ≥ 0") data.func(
                    λ;
                    derivative = -1,
                )
            end

            λ = 0.01
            @test_throws ArgumentError data.func(λ)

            @testset "Vectorization" begin
                λs = [0.4, 0.8, 1.2]
                @test data.func.(λs) isa Vector{Float64}
                @test issorted(data.func.(λs), rev = true)
            end
        end
    end
end
