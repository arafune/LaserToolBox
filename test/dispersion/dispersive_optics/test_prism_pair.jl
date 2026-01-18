using Test
using LaserToolBox
const ri = LaserToolBox.n

const DO = LaserToolBox.Dispersion.DispersiveOptics
# small tolerance for floating comparisons
const TOL = 1e-8

# ------ Short & simplest test for PrismPair ------

# SF11
# @ 800nm
n = ri.sf11(0.8)
brewster_angle = brewster_angle_deg(n)
ideal_apex = ideal_apex_deg(n)

@test brewster_angle_deg(ri.sf11, 0.8) ≈ brewster_angle atol=TOL
pr_brewster = PrismPair(
    wavelength = 800.0,
    incident_angle = brewster_angle,
    apex_angle = ideal_apex,
    insertion = (5.0, 10.0),
    material = ri.sf11,
)

@testset "Ideal Prism Pair (Brewster angle incidence)" begin
    prism = pr_brewster
    @test rad2deg.(DO.theta_1(prism)) ≈ rad2deg.(DO.theta_2(prism)) atol=TOL

    @test sin(deg2rad(brewster_angle)) ≈ n / sqrt(1+n^2) atol=TOL
    @test cos(deg2rad(brewster_angle)) ≈ 1.0 / sqrt(1+n^2) atol=TOL

    @test DO.dθ1_dΩ(prism) ≈ -(1.0/n^2) * DO.dn_dΩ(prism) atol=TOL
    @test DO.dθ3_dΩ(prism) ≈ 2DO.dn_dΩ(prism) atol=TOL
    @test DO.dn_dΩ(prism)[1] > 0

    l1 = prism.insertion[1] * cos(DO.theta_1(prism)[1]) / cos(DO.theta_2(prism)[1])
    @test DO.lg(prism)[1] ≈ (l1 + prism.insertion[2]) * sin(prism.apex_angle / 2) * 2
end


@testset "PrismPair constructor and wavelength normalization" begin
    pr1 = PrismPair(wavelength = 800.0)          # 800 nm -> 0.8 µm internally
    @test isa(pr1.wavelength, AbstractVector{Float64})
    @test length(pr1.wavelength) == 1

    pr2 = PrismPair(wavelength = [750.0, 800.0]) # nm -> µm vectorized
    @test isa(pr2.wavelength, AbstractVector{Float64})
    @test length(pr2.wavelength) == 2
end

@testset "Theta, lg and derivative helpers shapes" begin
    pr = PrismPair(wavelength = [800.0, 820.0])
    θ1 = DO.theta_1(pr)
    θ2 = DO.theta_2(pr)
    Lg = DO.lg(pr)
    @test length(θ1) == length(pr.wavelength)
    @test length(θ2) == length(pr.wavelength)
    @test length(Lg) == length(pr.wavelength)

    # Also check derivatives produce same-length vectors
    dn = DO.dn_dΩ(pr)
    dθ1 = DO.dθ1_dΩ(pr)
    dθ3 = DO.dθ3_dΩ(pr)
    @test length(dn) == length(pr.wavelength)
    @test length(dθ1) == length(pr.wavelength)
    @test length(dθ3) == length(pr.wavelength)
end

@testset "GDD: positive, negative and total" begin
    pr_single = PrismPair(wavelength = 800.0)
    gp_single = DO.gdd_positive(pr_single)
    gn_single = DO.gdd_negative(pr_single)
    total_single = gdd(pr_single)

    # For single-wavelength input, internal results should be vectors of length 1
    @test length(gp_single) == 1
    @test length(gn_single) == 1
    # gdd wrapper should return a scalar for single-wavelength input
    @test isa(total_single, Number)
    @test isapprox(total_single, gp_single[1] + gn_single[1]; atol = TOL, rtol = TOL)

    pr_multi = PrismPair(wavelength = [750.0, 800.0, 850.0])
    gp = DO.gdd_positive(pr_multi)
    gn = DO.gdd_negative(pr_multi)
    total = gdd(pr_multi)

    @test length(gp) == length(pr_multi.wavelength)
    @test length(gn) == length(pr_multi.wavelength)
    @test isa(total, AbstractVector)
    @test length(total) == length(pr_multi.wavelength)
    @test all(isapprox.(total, gp .+ gn; atol = TOL, rtol = TOL))
end

@testset "Brewster angle and ideal apex angle" begin
    # Use material function from an instance to avoid depending on module paths
    pr = PrismPair(wavelength = 800.0)
    λ = pr.wavelength[1]  # µm
    mat = pr.material

    # Compute via API
    b = brewster_angle_deg(mat, λ)
    ia = ideal_apex_deg(mat, λ)

    @test isfinite(b) && b > 0
    @test isfinite(ia) && ia > 0

    # Compare to direct formula using refractive index
    n = mat(λ)  # refractive index at λ (should be scalar)
    @test isapprox(b, rad2deg(atan(n)); atol = TOL, rtol = TOL)
end
