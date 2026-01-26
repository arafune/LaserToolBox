```julia
"""
    WavePlate(Γ::Real, θ::Real=0.0)

General wave plate with phase retardation Γ and fast axis at angle θ (in radians).

# Arguments
- `Γ::Real`: Phase retardation (in radians)
  - `Γ = π/2` for quarter-wave plate
  - `Γ = π` for half-wave plate
- `θ::Real`: Angle of the fast axis in radians (default: 0.0)

# Examples
```julia
# Quarter-wave plate with fast axis horizontal
Q = WavePlate(π/2, 0)

# Half-wave plate with fast axis at 22.5°
H = WavePlate(π, π/8)

# Arbitrary retardation
W = WavePlate(π/3, π/4)
```

# Notes
The Jones matrix for a general wave plate is:
```
M(Γ, θ) = R(-θ) · [1  0; 0  exp(iΓ)] · R(θ)
```
where R(θ) is a rotation matrix.
"""
function WavePlate(Γ::Real, θ::Real=0.0)
    c = cos(θ)
    s = sin(θ)
    
    # Phase factors
    eiΓ = exp(1im * Γ)
    
    # Jones matrix: R(-θ) * [1 0; 0 exp(iΓ)] * R(θ)
    M = [c^2 + eiΓ*s^2           c*s*(1 - eiΓ);
         c*s*(1 - eiΓ)           s^2 + eiΓ*c^2]
    
    return JonesMatrix(M)
end

"""
    QuarterWavePlate(θ::Real=0.0)

Quarter-wave plate (λ/4 plate) with fast axis at angle θ (in radians).

Equivalent to `WavePlate(π/2, θ)`.

# Arguments
- `θ::Real`: Angle of the fast axis in radians (default: 0.0)

# Examples
```julia
# QWP with fast axis horizontal
Q = QuarterWavePlate(0)

# Convert horizontal to right circular
j_out = Q * horizontal()
```
"""
QuarterWavePlate(θ::Real=0.0) = WavePlate(π/2, θ)

"""
    HalfWavePlate(θ::Real=0.0)

Half-wave plate (λ/2 plate) with fast axis at angle θ (in radians).

Equivalent to `WavePlate(π, θ)`.

# Arguments
- `θ::Real`: Angle of the fast axis in radians (default: 0.0)

# Examples
```julia
# HWP with fast axis at 22.5° rotates horizontal to vertical
H = HalfWavePlate(π/8)
j_out = H * horizontal()
```
"""
HalfWavePlate(θ::Real=0.0) = WavePlate(π, θ)
