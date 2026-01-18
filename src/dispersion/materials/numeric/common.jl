const DISPERSION_COMMON = """
# Arguments
 - `λ::Real`: Wavelength in μm

# Keyword Arguments
 - `derivative::Int=0`: Order of derivative to compute
"""

const DISPERSION_DOC = """ 
$DISPERSION_COMMON

# Returns
- `n::Float64`: Refractive index at wavelength `λ`

"""

const DISPERSION_BIAXIAL_DOC = """
$DISPERSION_COMMON

# Returns
  
- (n_e::Float64, n_o::Float64): NamedTuple where
    - n_e: extraordinary refractive index at wavelength λ
    - n_o: ordinary refractive index at wavelength λ

# Notes
- The extraordinary and ordinary refractive indices are defined with respect to the optical axis of the crystal.
- _e and _o functions are also provided for convenience.
"""
