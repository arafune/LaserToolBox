using ForwardDiff

"""
    beta_n(refractive_index_model, λ, n::Int)

Calculate the n-th order dispersion coefficient (βₙ) recursively.
- n=1: Group delay per unit length (1/v_g)
- n=2: Group Velocity Dispersion (GVD)
- n=3: Third-Order Dispersion (TOD)
- n=4: Fourth-Order Dispersion (FOD)

This function leverages ForwardDiff's ability to nest dual numbers, 
effectively calculating higher-order derivatives of the refractive index 
provided by the model's `derivative` argument.
"""
function beta_n(refractive_index_model, λ, derivative::Int = 0)
    c = 0.299792458  # µm/fs
    if n == 1
        # Base case: β₁ = dβ/dω = (n - λ * dn/dλ) / c
        # The model function is expected to support: model(λ, derivative=k)
        n_val = refractive_index_model(λ, derivative = 0)
        dn_dλ = refractive_index_model(λ, derivative = 1)
        return (n_val - λ * dn_dλ) / c
    else
        # Recursive step: βₙ = dβₙ₋₁/dω = (-λ² / 2πc) * dβₙ₋₁/dλ
        # We define βₙ₋₁ as a function of λ and differentiate it.
        prev_order(l) = beta_n(model, l, n-1)

        prefactor = -λ^2 / (2π * c)
        return prefactor * ForwardDiff.derivative(prev_order, λ)
    end
end

