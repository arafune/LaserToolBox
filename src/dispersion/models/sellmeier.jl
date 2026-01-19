module Sellmeier

using LaserToolBox.Derivatives: nth_derivative

"""
  sellmeier_expr(λ, A, B, C)

Expression for the Sellmeier equation.
"""
sellmeier_expr(λ, A, B, C) = begin
    λ2 = λ^2
    sqrt(
        one(λ) +
        A +
        B[1] * λ2 / (λ2 - C[1]) +
        B[2] * λ2 / (λ2 - C[2]) +
        B[3] * λ2 / (λ2 - C[3]),
    )
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
    B::NTuple{3,TCoeff},
    C::NTuple{3,TCoeff};
    derivative::Int = 0,
) where {T<:Real, TCoeff<:Real}
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

end
