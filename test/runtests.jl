using LaserToolBox
using Test

@testset "LaserToolBox.jl" begin
    @testset "Dispersion" begin
        @testset "Materials" begin
            @testset "Numeric" begin
                include("dispersion/materials/numeric/test_air.jl")
                include("dispersion/materials/numeric/test_glass.jl")
                include("dispersion/materials/numeric/test_biaxial.jl")
                include("dispersion/materials/numeric/test_bbo.jl")
                #
                include("dispersion/dispersive_optics/test_prism_pair.jl")
                include("dispersion/dispersive_optics/test_prism_pair_internal.jl")
                #
                include("dispersion/test_orders.jl")
                #
                #                include("optics/test_abcd.jl")
            end
        end
    end
end
