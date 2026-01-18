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

# @pytest ではなく @testset を使用し、verbose=true で詳細を表示
@testset verbose = true "Biaxial Material Tests" begin
    for data in biaxial_material_test_data
        # 各材料ごとのテストセット
        @testset "$(data.label) (Numeric)" begin
            λ = data.λ
            n_calc = data.func(λ)

            # expected[1] が n_e, expected[2] が n_o と仮定
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

