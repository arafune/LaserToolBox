#dispersion/dispersion.jl
module Dispersion
include("orders.jl")

export beta_n

"""
  `gvd`(model, λ; unit = :mm)


Group Velocity Dispersion (GVD) in [fs²/μm]
"""
function gvd(model, λ; unit = :mm)
    beta_n(model, λ; order = 2, unit)
end

"""

  `tod`(model, λ; unit = :mm)

Third-Order Dispersion (TOD) in [fs³/μm]
"""
function tod(model, λ; unit = :mm)
    beta_n(model, λ; order = 3, unit)
end

export gvd, tod

include("materials/materials.jl")
using .Materials
using ..Models

include("dispersive_optics/dispersive_optics.jl")
using .DispersiveOptics: PrismPair, gdd, brewster_angle_deg, ideal_apex_deg

export DispersiveOptics
export PrismPair, gdd, brewster_angle_deg, ideal_apex_deg

import ..RefractiveIndex
const n = RefractiveIndex
end
