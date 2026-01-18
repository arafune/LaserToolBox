# dispersion/materials/numeric/glass.jl

using ..MaterialConstants: BK7, BK7_SOURCE
using ..MaterialConstants: FusedSilica, FUSED_SILICA_SOURCE
using ..MaterialConstants: CaF2, CaF2_SOURCE
using ..MaterialConstants: SF10, SF10_SOURCE
using ..MaterialConstants: SF11, SF11_SOURCE

using ...Models.Sellmeier: sellmeier

"""
  bk7(λ; derivative=0)

Dispersion of BK7.

$BK7_SOURCE

$DISPERSION_DOC
"""
function bk7(λ; derivative::Int = 0)
    if λ < BK7.range[1] || λ>BK7.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(BK7.range[1]) - $(BK7.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, BK7.A, BK7.B, BK7.C; derivative = derivative)
end

"""
  caf2(λ; derivative=0)

Dispersion of CaF2 (0.15-12 μm).

$CaF2_SOURCE

$DISPERSION_DOC
"""
function caf2(λ; derivative::Int = 0)
    if λ < CaF2.range[1] || λ>CaF2.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(CaF2.range[1]) - $(CaF2.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, CaF2.A, CaF2.B, CaF2.C; derivative = derivative)
end

"""
  fused_silica(λ; derivative=0)

Dispersion of Fused Silica.

$FUSED_SILICA_SOURCE

$DISPERSION_DOC
"""
function fused_silica(λ; derivative::Int = 0)
    if λ < FusedSilica.range[1] || λ>FusedSilica.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(FusedSilica.range[1]) - $(FusedSilica.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(
        λ,
        FusedSilica.A,
        FusedSilica.B,
        FusedSilica.C;
        derivative = derivative,
    )
end


"""
  sf10(λ; derivative=0)

Dispersion of SF10.

$SF10_SOURCE

$DISPERSION_DOC
"""
function sf10(λ; derivative::Int = 0)
    if λ < SF10.range[1] || λ>SF10.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(SF10.range[1]) - $(SF10.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, SF10.A, SF10.B, SF10.C; derivative = derivative)
end

"""
  sf11(λ; derivative=0)

Dispersion of SF11.

$SF11_SOURCE

$DISPERSION_DOC
"""
function sf11(λ; derivative::Int = 0)
    if λ < SF11.range[1] || λ>SF11.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(SF11.range[1]) - $(SF11.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, SF11.A, SF11.B, SF11.C; derivative = derivative)
end

