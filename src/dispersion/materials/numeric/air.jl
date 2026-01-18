# dispersion/materials/numeric/air.jl

using ...Models.AirDispersion: air_dispersion
using ..MaterialConstants: AIR_SOURCE, Air

"""
  air(λ; derivative=0)

Dispersion of Air.

$AIR_SOURCE

$DISPERSION_DOC
"""
function air(λ; derivative::Int = 0)
    if λ < Air.range[1] || λ>Air.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(Air.range[1]) - $(Air.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return air_dispersion(λ, Air.B, Air.C; derivative = derivative)
end
