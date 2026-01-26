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
export Models

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

include("optics/abcd.jl")

export Medium, FreeSpace, ThinLens, ThickLens, Interface, PlanoConvexLens
export transfer_matrix, effective_focal_length, back_focal_length


include("polarization/jones.jl")

export JonesVector, JonesMatrix
export LinearPolarizer, QuarterWavePlate, HalfWavePlate, Rotator, WavePlate
export horizontal, vertical, diagonal, antidiagonal, right_circular, left_circular
export intensity, degree_of_polarization, ellipse_parameters


end
