#dispersion/dispersion.jl
module Dispersion
include("materials/materials.jl")
include("dispersive_optics/dispersive_optics.jl")
include("orders.jl")

"""Group Velocity Dispersion (GVD) in [fs²/μm]"""
const gvd(model, λ) = beta_n(model, λ, 2)

"""Third-Order Dispersion (TOD) in [fs³/μm]"""
const tod(model, λ) = beta_n(model, λ, 3)

using .Materials
using ..Models

using .DispersiveOptics: PrismPair, gdd, brewster_angle_deg, ideal_apex_deg

export DispersiveOptics
export PrismPair, gdd, brewster_angle_deg, ideal_apex_deg

export beta_n, gvd, tod

import ..RefractiveIndex
const n = RefractiveIndex
end
