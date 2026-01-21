using Test
using LaserToolBox


"""
Edumund #11-720 refractive index tests
tc =1.6mm, R = 114.75 mm
@587.6 mm f= 249.19mm
"""
@testset "Edumund #11-720" begin
    edmund_11_720 = PlanoConvexLens(n_lens = n.fused_silica(0.5876), d = 1.6, R = 114.75)
    @test back_focal_length(edmund_11_720) ≈ 249.19 atol=0.1

end

