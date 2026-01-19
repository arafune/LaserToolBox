## LaserToolBox.jl
#
# Main module file for the LaserToolBox package.
# This package provides tools for numerical differentiation, dispersion modeling,
# material constants, and refractive index calculations relevant to laser spectroscopy.

module LaserToolBox

include("numerics/derivatives.jl")
using .Derivatives: nth_derivative

include("dispersion/models/models.jl")
using .Models

include("dispersion/materials/constants/constants.jl")
using .MaterialConstants

include("dispersion/materials/numeric/numeric.jl")
import .RefractiveIndex
export RefractiveIndex
const n = RefractiveIndex
export n

include("dispersion/dispersion.jl")

using .Dispersion

export Dispersion
export PrismPair, gdd, brewster_angle_deg, ideal_apex_deg
export beta_n, gvd, tod

end
