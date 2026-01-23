module Sellmeier
using Symbolics

using LaserToolBox.Derivatives: nth_derivative

"""
  sellmeier_expr(λ, A, B, C)

Expression for the Sellmeier equation.
"""
sellmeier_expr(λ, A, B, C) = begin
    λ2 = λ^2
    sqrt(one(λ) + A + sum(B[i] * λ2 / (λ2 - C[i]) for i in eachindex(B, C)))
end


"""
  sellmeier(λ, A, B, C; derivative=0)

Calculate the refractive index using the Sellmeier equation.

# Arguments

- λ: Wavelength in micrometers (μm)
- A: Coefficient A
- B: Coefficients B_i
- C: Coefficients C_i
- derivative: Order of derivative to compute (default is 0)

# Returns
- Refractive index n at wavelength λ
"""
function sellmeier(
    λ::T,
    A::TCoeff,
    B::Tuple{Vararg{TCoeff}},
    C::Tuple{Vararg{TCoeff}};
    derivative::Int = 0,
) where {T<:Real,TCoeff<:Real}
    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    if derivative == 0
        return sellmeier_expr(λ, A, B, C)
    else
        f(x) = sellmeier_expr(x, A, B, C)
        return nth_derivative(f, λ, derivative)
    end
end

"""
sellmeier_sym(derivative=0)
"""
function sellmeier_sym(derivative::Int = 0)
    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    @variables λ A B[1:3] C[1:3]

    n = sellmeier_expr(λ, A, B, C)

    return derivative == 0 ? n : (Differential(λ)^derivative) * n
end

export sellmeier_sym
end
