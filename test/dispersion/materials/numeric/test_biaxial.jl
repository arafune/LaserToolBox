using Test
using LaserToolBox

# const n = LaserToolBox.n 

biaxial_material_test_data = [
    (func = n.mgf2, expected = (1.3896, 1.3777), atol = 1e-3, λ = 0.5876, label = "MgF2")
    (
        func = n.calcite,
        expected = (1.480, 1.642),
        atol = 1e-3,
        λ = 1.064,
        label = "Calcite",
    )
    (
        func = n.quartz,
        expected = (1.5472301086112594, 1.5383355123424691),
        atol = 1e-13,
        λ = 0.8,
        label = "Quartz",
    )
]

@testset verbose = true "Biaxial Material Tests" begin
    for data in biaxial_material_test_data

        @testset "$(data.label) (Numeric)" begin
            λ = data.λ
            n_calc = data.func(λ)

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



biaxial_material_e_axis_test_data = [
    (func = n.mgf2_e, expected = 1.3896, atol = 1e-3, λ = 0.5876, label = "MgF2")
    (func = n.calcite_e, expected = 1.480, atol = 1e-3, λ = 1.064, label = "Calcite")
    (
        func = n.quartz_e,
        expected = 1.5472301086112594,
        atol = 1e-13,
        λ = 0.8,
        label = "Quartz",
    )
]


@testset "BiAxial material Refractive Indices along e-axis (Numeric)" begin
    for data in biaxial_material_e_axis_test_data
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


biaxial_material_o_axis_test_data = [
    (func = n.mgf2_o, expected = 1.3777, atol = 1e-3, λ = 0.5876, label = "MgF2")
    (func = n.calcite_o, expected = 1.642, atol = 1e-3, λ = 1.064, label = "Calcite")
    (
        func = n.quartz_o,
        expected = 1.5383355123424691,
        atol = 1e-13,
        λ = 0.8,
        label = "Quartz",
    )
]


@testset "BiAxial material Refractive Indices along o-axis (Numeric)" begin
    for data in biaxial_material_o_axis_test_data
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
