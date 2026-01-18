using Test
using LaserToolBox
const ri = LaserToolBox.n

const DO = LaserToolBox.Dispersion.DispersiveOptics
#using LaserToolBox.Dispersion: theta_1, theta_2, lg,
#                               dn_dΩ, dθ1_dΩ, dθ3_dΩ,
# small tolerance for floating comparisons
const TOL = 1e-8


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

# Theta / GDD / Brewster tests using the DispersiveOptics namespace
pr = PrismPair(wavelength = [800.0, 820.0])
θ1 = DO.theta_1(pr)
θ2 = DO.theta_2(pr)
Lg = DO.lg(pr)

@testset "Theta, lg and derivative helpers shapes" begin
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
