module ModelSymbolic

using Symbolics

"""
air_dispersion_sym(derivative=0)
"""
function air_dispersion_sym(derivative::Int = 0)
    if derivative < 0
        throw(ArgumentError("derivative must be ≥ 0"))
    end

    @variables λ B[1:2] C[1:2]

    n = air_dispersion_expr(λ, B, C)

    return derivative == 0 ? n : Differential(λ)^derivative(n)
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

    return derivative == 0 ? n : Differential(λ)^derivative(n)
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

    return derivative == 0 ? n : Differential(λ)^derivative(n)
end

end
