# dispersion/materials/numeric/biaxial.jl
using ..MaterialConstants: MgF2, MgF2_SOURCE
using ..MaterialConstants: Calcite, CALCITE_SOURCE
using ..MaterialConstants: Quartz, QUARTZ_SOURCE

using ...Models.Sellmeier: sellmeier

"""
calcite_e(λ; derivative=0)

Dispersion of the extraordinary index of Calcite (CaCO3).

$CALCITE_SOURCE

$DISPERSION_DOC
"""
function calcite_e(λ; derivative::Int = 0)
    if λ < Calcite.range[1] || λ>Calcite.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(Calcite.range[1]) - $(Calcite.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, Calcite.e.A, Calcite.e.B, Calcite.e.C; derivative = derivative)
end

"""
calcite_o(λ; derivative=0)

Dispersion of the ordinary index of Calcite (CaCO3).

$CALCITE_SOURCE

$DISPERSION_DOC
"""
function calcite_o(λ; derivative::Int = 0)
    if λ < Calcite.range[1] || λ>Calcite.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(Calcite.range[1]) - $(Calcite.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, Calcite.o.A, Calcite.o.B, Calcite.o.C; derivative = derivative)
end



"""
  calcite(λ; derivative=0)

Dispersion of Calcite (CaCO3).

$CALCITE_SOURCE

$DISPERSION_BIAXIAL_DOC
"""
function calcite(λ; derivative::Int = 0)
    n_e = calcite_e(λ, ; derivative = derivative)
    n_o = calcite_o(λ, ; derivative = derivative)
    return (n_e = n_e, n_o = n_o)
end

"""
  mgf2_e(λ; derivative=0)

Dispersion of the extraordinary index of MgF₂.

$MgF2_SOURCE

$DISPERSION_DOC
"""
function mgf2_e(λ; derivative::Int = 0)
    if λ < MgF2.range[1] || λ>MgF2.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(MgF2.range[1]) - $(MgF2.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, MgF2.e.A, MgF2.e.B, MgF2.e.C; derivative = derivative)
end

"""
  mgf2_o(λ; derivative=0)
Dispersion of the ordinary index of MgF₂.

$MgF2_SOURCE

$DISPERSION_DOC
"""
function mgf2_o(λ; derivative::Int = 0)
    if λ < MgF2.range[1] || λ>MgF2.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(MgF2.range[1]) - $(MgF2.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, MgF2.o.A, MgF2.o.B, MgF2.o.C; derivative = derivative)
end

"""
  mgf2(λ; derivative=0)

Dispersion of MgF₂ (biaxial crystal). 

$MgF2_SOURCE

$DISPERSION_BIAXIAL_DOC
"""
function mgf2(λ; derivative::Int = 0)
    n_e = mgf2_e(λ, ; derivative = derivative)
    n_o = mgf2_o(λ, ; derivative = derivative)
    return (n_e = n_e, n_o = n_o)
end

"""
  quartz_e(λ; derivative=0)
Dispersion of the extraordinary index of crystal quartz.

$QUARTZ_SOURCE

$DISPERSION_DOC 
"""
function quartz_e(λ; derivative::Int = 0)
    if λ < Quartz.range[1] || λ>Quartz.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(Quartz.range[1]) - $(Quartz.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, Quartz.e.A, Quartz.e.B, Quartz.e.C; derivative = derivative)
end

"""
  quartz_o(λ; derivative=0)

Dispersion of the ordinary index of crystal quartz. 
$QUARTZ_SOURCE

$DISPERSION_DOC
"""
function quartz_o(λ; derivative::Int = 0)
    if λ < Quartz.range[1] || λ>Quartz.range[2]
        msg = "Wavelength $λ µm is out of the valid range for dispersion model ($(Quartz.range[1]) - $(Quartz.range[2]) µm)."
        throw(ArgumentError(msg))
    end
    return sellmeier(λ, Quartz.o.A, Quartz.o.B, Quartz.o.C; derivative = derivative)
end



"""
  quartz(λ; derivative=0)
  
Dispersion of crystal quartz.

$QUARTZ_SOURCE
 
$DISPERSION_BIAXIAL_DOC
"""
function quartz(λ; derivative::Int = 0)
    n_e = quartz_e(λ, ; derivative = derivative)
    n_o = quartz_o(λ, ; derivative = derivative)
    return (n_e = n_e, n_o = n_o)
end

