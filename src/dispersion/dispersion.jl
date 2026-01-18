#dispersion/dispersion.jl
module Dispersion
include("materials/materials.jl")
include("dispersive_optics/dispersive_optics.jl")

using .Materials
using ..Models

using .DispersiveOptics: PrismPair, gdd, brewster_angle_deg, ideal_apex_deg

export DispersiveOptics
export PrismPair, gdd, brewster_angle_deg, ideal_apex_deg
end
