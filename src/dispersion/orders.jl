using ForwardDiff

using .Materials.RefractiveIndex
"""
    beta_n(refractive_index_model, λ; order::Int = 1, unit::Symbol = :mm)

Calculate the n-th order dispersion coefficient (βₙ) recursively for a given refractive index model and wavelength.

# Arguments
- `refractive_index_model`: Function that returns the refractive index and its derivatives. Must accept keyword argument `derivative`.
- `λ`: Wavelength in micrometers (μm).
- `order::Int=1`: Order of the dispersion coefficient (βₙ). 
    - 1: Group delay per unit length (1/v_g)
    - 2: Group Velocity Dispersion (GVD)
    - 3: Third-Order Dispersion (TOD)
    - 4: Fourth-Order Dispersion (FOD)
- `unit::Symbol=:mm`: Output length unit. `:mm` for per millimeter (default), `:μm` for per micrometer.

# Returns
- The n-th order dispersion coefficient βₙ in units of fsⁿ/(unit).

# Notes
- This function uses ForwardDiff to compute higher-order derivatives with respect to wavelength.
- The result is automatically scaled to the requested length unit.

# Example
```julia
beta_n(n.sf11, 0.8; order=2, unit=:mm)  # GVD in fs^2/mm
```
"""
function beta_n(refractive_index_model, λ; order::Int = 1, unit::Symbol = :mm)
    c = 0.299792458  # µm/fs
    scale = unit == :mm ? 1000^(order-1) : 1.0  # Scale factor for length units
    if order== 1
        # Base case: β₁ = dβ/dω = (n - λ * dn/dλ) / c
        # The model function is expected to support: model(λ, derivative=k)
        n_val = refractive_index_model(λ; derivative = 0)
        dn_dλ = refractive_index_model(λ; derivative = 1)
        return (n_val - λ * dn_dλ) / c * scale
    else
        # Recursive step: βₙ = dβₙ₋₁/dω = (-λ² / 2πc) * dβₙ₋₁/dλ
        # We define βₙ₋₁ as a function of λ and differentiate it.
        prev_order(l) = beta_n(refractive_index_model, l; order=order-1)

        prefactor = -λ^2 / (2π * c)
        return prefactor * ForwardDiff.derivative(prev_order, λ) * scale
    end
end

