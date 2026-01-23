using LaserToolBox
using Coverage
using Test
## if Base.JLOptions().code_coverage == 0
##     println("Restarting with --code-coverage for coverage measurement...")
##     run(
##         `$(Base.julia_cmd()) --project=$(Base.active_project()) --code-coverage test/runtests.jl`,
##     )
##     exit()
## end

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
                #include("dispersion/models/test_model_symbolic.jl")
                #
                include("optics/test_abcd.jl")
            end
        end
    end
end


results = process_folder()
LCOV.writefile("lcov.info", results)
run(`genhtml -o coverage_html lcov.info`)
