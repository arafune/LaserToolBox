using Test
using Symbolics: Num, ComposedFunction
using LaserToolBox

@testset "sellmeler_sym" begin
    @testset "Invalid Derivative" begin
        @test_throws ArgumentError("derivative must be ≥ 0") LaserToolBox.Dispersion.Models.Sellmeier.sellmeier_sym(
            -1,
        )
    end

    @testset "Zeroth Derivative" begin
        refractive_index = LaserToolBox.Dispersion.Models.Sellmeier.sellmeier_sym(0)
        @test typeof(refractive_index) == Num
    end

    @testset "First Derivative" begin
        dn = LaserToolBox.Dispersion.Models.Sellmeier.sellmeier_sym(1)
        @test dn isa ComposedFunction || dn isa Num
    end

    @testset "Second Derivative" begin
        d2n = LaserToolBox.Dispersion.Models.Sellmeier.sellmeier_sym(2)
        @test d2n isa ComposedFunction || dwn isa Num
    end
end


@testset "bbo_sellmeier_sym" begin
    @testset "Invalid Derivative" begin
        @test_throws ArgumentError("derivative must be ≥ 0") LaserToolBox.Dispersion.Models.BBOSellmeier.bbo_sellmeier_sym(
            -1,
        )
    end

    @testset "Zeroth Derivative" begin
        refractive_index = LaserToolBox.Dispersion.Models.BBOSellmeier.bbo_sellmeier_sym(0)
        @test typeof(refractive_index) == Num
    end

    @testset "First Derivative" begin
        dn = LaserToolBox.Dispersion.Models.BBOSellmeier.bbo_sellmeier_sym(1)
        @test dn isa ComposedFunction || dn isa Num
    end

    @testset "Second Derivative" begin
        d2n = LaserToolBox.Dispersion.Models.BBOSellmeier.bbo_sellmeier_sym(2)
        @test d2n isa ComposedFunction || dwn isa Num
    end
end

@testset "air_dispersion_sym" begin
    @testset "Invalid Derivative" begin
        @test_throws ArgumentError("derivative must be ≥ 0") LaserToolBox.Dispersion.Models.AirDispersion.air_dispersion_sym(
            -1,
        )
    end

    @testset "Zeroth Derivative" begin
        refractive_index =
            LaserToolBox.Dispersion.Models.AirDispersion.air_dispersion_sym(0)
        @test typeof(refractive_index) == Num
    end

    @testset "First Derivative" begin
        dn = LaserToolBox.Dispersion.Models.AirDispersion.air_dispersion_sym(1)
        @test dn isa ComposedFunction || dn isa Num
    end

    @testset "Second Derivative" begin
        d2n = LaserToolBox.Dispersion.Models.AirDispersion.air_dispersion_sym(2)
        @test d2n isa ComposedFunction || dwn isa Num
    end
end
