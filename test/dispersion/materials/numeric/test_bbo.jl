using Test
using LaserToolBox

# const n = LaserToolBox.n 

biaxial_material_test_data = [
    (
        func = n.alpha_bbo,
        expected = (1.534, 1.673),
        atol = 1e-3,
        λ = 0.5876,
        label = "Alpha BBO d-line",
    )
    (
        func = n.alpha_bbo,
        expected = (1.5284, 1.6639),
        atol = 1e-3,
        λ = 0.8,
        label = "Alpha BBO 800nm",
    )
    (
        func = n.alpha_bbo,
        expected = (1.5500, 1.6962),
        atol = 1e-3,
        λ = 0.4,
        label = "Alpha BBO 400nm",
    )
    (
        func = n.beta_bbo,
        expected = (1.5462, 1.6614),
        atol = 1e-3,
        λ = 0.8,
        label = "beta BBO 800nm",
    )
]

# 
@testset verbose = true "Biaxial Material Tests" begin
    for data in biaxial_material_test_data
        #
        @testset "$(data.label) (Numeric)" begin
            λ = data.λ
            n_calc = data.func(λ)

            # Assumption expected[1] is n_e, expected[2] is n_o
            @test n_calc.n_e ≈ data.expected[1] atol = data.atol
            @test n_calc.n_o ≈ data.expected[2] atol = data.atol

            @testset "Derivative Consistency" begin
                dn_e = data.func(λ; derivative = 1).n_e
                dn_o = data.func(λ; derivative = 1).n_o
                @test dn_e < 0
                @test dn_o < 0

                d2n_e = data.func(λ; derivative = 2).n_e
                d2n_o = data.func(λ; derivative = 2).n_o
                @test d2n_e > 0
                @test d2n_o > 0

                # Test Error message
                @test_throws ArgumentError("derivative must be ≥ 0") data.func(
                    λ;
                    derivative = -1,
                ).n_e
                @test_throws ArgumentError("derivative must be ≥ 0") data.func(
                    λ;
                    derivative = -1,
                ).n_o

            end

            λ = 0.01
            @test_throws ArgumentError data.func(λ)

            @testset "Vectorization" begin
                λs = [0.4, 0.8, 1.2]
                n_results = data.func.(λs)
                n_e_vals = getproperty.(n_results, :n_e)
                n_o_vals = getproperty.(n_results, :n_o)

                @test n_e_vals isa Vector{Float64}
                @test n_o_vals isa Vector{Float64}
                @test issorted(n_e_vals, rev = true)
                @test issorted(n_o_vals, rev = true)
            end
        end
    end
end


bbo_e_axis_test_data = [
    (
        func = n.beta_bbo_e,
        expected = 1.5462,
        atol = 1e-3,
        λ = 0.8,
        label = "beta BBO 800nm (e)",
    )
    (
        func = n.alpha_bbo_e,
        expected = 1.5500,
        atol = 1e-3,
        λ = 0.4,
        label = "Alpha BBO 400nm (e)",
    )
]

@testset "BBO Refractive Indices  along e-axis (Numeric)" begin
    for data in bbo_e_axis_test_data
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
                λs = [0.4, 0.8, 1.1]
                @test data.func.(λs) isa Vector{Float64}
                @test issorted(data.func.(λs), rev = true)
            end
        end
    end
end
