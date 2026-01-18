# numeric/numeric.jl
# This file collects all numeric computation functions for materials
# and exposes them via a single module.

module RefractiveIndex

include("common.jl")

using LaserToolBox.Derivatives: nth_derivative

using ..MaterialConstants

# Include numeric implementations for each material
include("air.jl")
export air

include("bbo.jl")
export alpha_bbo, beta_bbo
export alpha_bbo_e, alpha_bbo_o
export beta_bbo_e, beta_bbo_o

include("biaxial.jl")
export calcite, mgf2, quartz
export calcite_e, calcite_o, mgf2_e, mgf2_o, quartz_e, quartz_o

include("glass.jl")
export bk7, caf2, fused_silica, sf10, sf11

end
