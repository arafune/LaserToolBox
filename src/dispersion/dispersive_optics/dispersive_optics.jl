# dispersion/dispersive_optics/dispersive_optics.jl

module DispersiveOptics
include("angles.jl")
include("prism_pair.jl")

export PrismPair, gdd, brewster_angle_deg, ideal_apex_deg
end

