# dispersion/materials/numeric/bbo.jl
using ..MaterialConstants: AlphaBBO, ALPHA_BBO_SOURCE
using ..MaterialConstants: BetaBBO, BETA_BBO_SOURCE

using ...Models.BBOSellmeier: bbo_sellmeier
using ...Models.Sellmeier: sellmeier

"""
  alpha_bbo_e(λ; derivative=0)

Dispersion of the extraordinary index of α-BBO.

$ALPHA_BBO_SOURCE

$DISPERSION_DOC
"""
function alpha_bbo_e(λ; derivative::Int = 0)
    if λ < AlphaBBO.range[1] || λ>AlphaBBO.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(AlphaBBO.range[1]) - $(AlphaBBO.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return bbo_sellmeier(λ, AlphaBBO.coeffs_e...; derivative = derivative)
end

"""
  alpha_bbo_o(λ; derivative=0)

Dispersion of the ordinary index of α-BBO.

$ALPHA_BBO_SOURCE

$DISPERSION_DOC
"""
function alpha_bbo_o(λ; derivative::Int = 0)
    if λ < AlphaBBO.range[1] || λ>AlphaBBO.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(AlphaBBO.range[1]) - $(AlphaBBO.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return bbo_sellmeier(λ, AlphaBBO.coeffs_o...; derivative = derivative)
end

"""
  alpha_bbo(λ; derivative=0)

Dispersion of α-BBO.

$ALPHA_BBO_SOURCE

$DISPERSION_BIAXIAL_DOC
"""
function alpha_bbo(λ; derivative::Int = 0)
    n_o = alpha_bbo_o(λ, derivative = derivative)
    n_e = alpha_bbo_e(λ, derivative = derivative)
    return (n_e = n_e, n_o = n_o)
end

"""
  beta_bbo_e(λ; derivative=0)

Dispersion of the extraordinary index of β-BBO.

$BETA_BBO_SOURCE

$DISPERSION_DOC
"""
function beta_bbo_e(λ; derivative::Int = 0)
    if λ < BetaBBO.range[1] || λ>BetaBBO.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(BetaBBO.range[1]) - $(BetaBBO.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, BetaBBO.e.A, BetaBBO.e.B, BetaBBO.e.C; derivative = derivative)
    #return bbo_sellmeier(λ, BetaBBO.coeffs_e...; derivative = derivative)
end

"""
  beta_bbo_o(λ; derivative=0)

Dispersion of the ordinary index of β-BBO.

$BETA_BBO_SOURCE

$DISPERSION_DOC
"""
function beta_bbo_o(λ; derivative::Int = 0)
    if λ < BetaBBO.range[1] || λ>BetaBBO.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(BetaBBO.range[1]) - $(BetaBBO.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, BetaBBO.o.A, BetaBBO.o.B, BetaBBO.o.C; derivative = derivative)
    #return bbo_sellmeier(λ, BetaBBO.coeffs_o...; derivative = derivative)
end

"""
  beta_bbo(λ; derivative=0)

Dispersion of β-BBO.

$BETA_BBO_SOURCE

$DISPERSION_BIAXIAL_DOC
"""
function beta_bbo(λ; derivative::Int = 0)
    n_o = beta_bbo_o(λ, derivative = derivative)
    n_e = beta_bbo_e(λ, derivative = derivative)
    return (n_e = n_e, n_o = n_o)
end
