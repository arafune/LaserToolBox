# dispersion/dispersive_optics/dispersive_optics.jl

module DispersiveOptics
include("angles.jl")
include("prism_pair.jl")

export PrismPair, gdd, gdd_positive, gdd_negative, lg, brewster_angle_deg, ideal_apex_deg
end
