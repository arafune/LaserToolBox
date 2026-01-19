# dispersion/models/ciddor.jl
module AirDispersion

using LaserToolBox.Derivatives: nth_derivative

"""
air_dispersion_expr(λ; derivative=0)

Ciddor-type formula, not Sellmeier.

"""
air_dispersion_expr(λ, B, C) = begin
    λm2 = λ^(-2)
    1 + B[1] / (C[1] - λm2) + B[2] / (C[2] - λm2)
end

"""
  air_dispersion(λ, B, C; derivative=0)

Calculate the refractive index using the Ciddor-type

# Arrguments

- λ: Wavelength in μm
- B: Coefficients B_i
- C: Coefficients C_i
- derivative: Order of derivative to compute (default is 0)

# Returns:
- Refractive index n at wavelength λ
"""
function air_dispersion(
    λ::T,
    B::NTuple{2,TCoeff},
    C::NTuple{2,TCoeff};
    derivative::Int = 0,
) where {T<:Real, TCoeff<:Real}
    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    if derivative == 0
        return air_dispersion_expr(λ, B, C)
    else
        f(x) = air_dispersion_expr(x, B, C)
        return nth_derivative(f, λ, derivative)
    end
end

end
