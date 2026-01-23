module BBOSellmeier
using Symbolics

using LaserToolBox.Derivatives: nth_derivative

"""
  bbo_sellmeier_expr(λ, A, B, C)

Expression for the Sellmeier equation.
"""
bbo_sellmeier_expr(λ, A, B, C, D) = begin
    λ2 = λ^2
    sqrt(A - D * λ2 + B / (λ2 - C))
end

"""
  bbo_sellmeier(λ, A, B, C; derivative=0)

Calculate the refractive index using the Sellmeier equation for BBO.

# Arguments

- λ: Wavelength in micrometers (μm)
- A: Coefficient A
- B: Coefficients h
- C: Coefficients C_i
- derivative: Order of derivative to compute (default is 0)

# Returns
- Refractive index n at wavelength λ
"""
function bbo_sellmeier(
    λ::T,
    A::TCoeff,
    B::TCoeff,
    C::TCoeff,
    D::TCoeff;
    derivative::Int = 0,
) where {T<:Real,TCoeff<:Real}
    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    if derivative == 0
        return bbo_sellmeier_expr(λ, A, B, C, D)
    else
        f(x) = bbo_sellmeier_expr(x, A, B, C, D)
        return nth_derivative(f, λ, derivative)
    end
end

"""
bbo_sellmeier_sym(derivative=0)
"""
function bbo_sellmeier_sym(derivative::Int = 0)

    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    @variables λ A B C D

    n = bbo_sellmeier_expr(λ, A, B, C, D)

    return derivative == 0 ? n : (Differential(λ)^derivative)*n
end

export bbo_sellmeier_sym
end
