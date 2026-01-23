#dispersion/dispersion.jl
module Dispersion
include("orders.jl")

export beta_n, gvd, tod

"""Group Velocity Dispersion (GVD) in [fs²/μm]"""
const gvd(model, λ; unit) = beta_n(model, λ; order = 2, unit = :mm)

"""Third-Order Dispersion (TOD) in [fs³/μm]"""
const tod(model, λ; unit) = beta_n(model, λ; order = 3, unit = :mm)

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
