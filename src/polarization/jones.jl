# polarization/jones.jl
using LinearAlgebra

export JonesVector, JonesMatrix
export LinearPolarizer, QuarterWavePlate, HalfWavePlate, Rotator
export horizontal, vertical, diagonal, antidiagonal, right_circular, left_circular
export intensity, degree_of_polarization, ellipse_parameters

"""
    JonesVector

Represents the polarization state of light as a complex 2-element vector.

# Fields
- `Ex::ComplexF64`: Electric field component in x-direction
- `Ey::ComplexF64`: Electric field component in y-direction

# Examples
```julia
# Horizontal polarization
h = JonesVector(1, 0)

# Right circular polarization
r = JonesVector(1, 1im) / sqrt(2)
```
"""
struct JonesVector
    Ex::ComplexF64
    Ey::ComplexF64
end

#Constructor from real numbers
JonesVector(Ex::Real, Ey::Real) = JonesVector(ComplexF64(Ex), ComplexF64(Ey))

# Vector representation
Base.Vector(j::JonesVector) = [j.Ex, j.Ey]
Base.Array(j::JonesVector) = Vector(j)

#Arithmetic operations
Base.:+(j1::JonesVector, j2::JonesVector) = JonesVector(j1.Ex + j2.Ex, j1.Ey + j2.Ey)
Base.:-(j1::JonesVector, j2::JonesVector) = JonesVector(j1.Ex - j2.Ex, j1.Ey - j2.Ey)
Base.:*(a::Number, j::JonesVector) = JonesVector(a * j.Ex, a * j.Ey)
Base.:*(j::JonesVector, a::Number) = a * j
Base.:/(j::JonesVector, a::Number) = JonesVector(j.Ex / a, j.Ey / a)

#Norm
LinearAlgebra.norm(j::JonesVector) = sqrt(abs2(j.Ex) + abs2(j.Ey))
LinearAlgebra.normalize(j::JonesVector) = j / norm(j)

#Intensity
"""
  intensity(j::JonesVector)

Calculate the intensity of the polarized light (proportional to |E|²).
"""
intensity(j::JonesVector) = abs2(j.Ex) + abs2(j.Ey)

"""
  JonesMatrix

Represents a polarization-transforming optical element as a 2×2 complex matrix.

# Fields
- M::Matrix{ComplexF64}: 2×2 Jones matrix

# Examples
```julia
# Linear polarizer at 0°
P = JonesMatrix([1 0; 0 0])

# Quarter-wave plate with fast axis at 0°
Q = JonesMatrix([1 0; 0 1im])
```
"""
struct JonesMatrix
    M::Matrix{ComplexF64}
end

#Constructor from real matrix
JonesMatrix(M::Matrix{<:Real}) = JonesMatrix(ComplexF64.(M))

#Matrix representation
Base.Matrix(j::JonesMatrix) = j.M
Base.Array(j::JonesMatrix) = Matrix(j)

# Matrix-vector multiplication (applying Jones matrix to Jones vector)
Base.:*(J::JonesMatrix, v::JonesVector) = begin
    result = J.M * Vector(v)
    JonesVector(result[1], result[2])
end

# Matrix-matrix multiplication (combining Jones matrices)
Base.:*(J1::JonesMatrix, J2::JonesMatrix) = JonesMatrix(J1.M * J2.M)

#Inverse
Base.inv(J::JonesMatrix) = JonesMatrix(inv(J.M))

#Determinant
LinearAlgebra.det(J::JonesMatrix) = det(J.M)

# ============================================================================
# Common polarization states
# ============================================================================

"""
  horizontal()

Horizontal linear polarization (x-direction).
"""
horizontal() = JonesVector(1, 0)

"""
  vertical()

Vertical linear polarization (y-direction).
"""
vertical() = JonesVector(0, 1)

"""
  diagonal()

Diagonal linear polarization (+45°).
"""
diagonal() = JonesVector(1, 1) / sqrt(2)

"""
  antidiagonal()

Anti-diagonal linear polarization (-45°).
"""
antidiagonal() = JonesVector(1, -1) / sqrt(2)

"""
  right_circular()

Right-handed circular polarization.
"""
right_circular() = JonesVector(1, 1im) / sqrt(2)

"""
  left_circular()

Left-handed circular polarization.
"""
left_circular() = JonesVector(1, -1im) / sqrt(2)

# ============================================================================
# Common optical elements
# ============================================================================

"""

  LinearPolarizer(θ::Real)

Linear polarizer with transmission axis at angle θ (in radians).

# Arguments
- θ::Real: Angle of the transmission axis in radians

# Examples
```julia
# Horizontal polarizer
P_h = LinearPolarizer(0)

# Vertical polarizer
P_v = LinearPolarizer(π/2)

# 45° polarizer
P_45 = LinearPolarizer(π/4)
```
"""
function LinearPolarizer(θ::Real)
    return Rotator(-θ) * JonesMatrix([1 0; 0 0]) * Rotator(θ)
end


"""
  Rotator(θ::Real)

Polarization rotator that rotates the polarization state by angle θ (in radians).

# Arguments
- θ::Real: Rotation angle in radians (positive is counterclockwise)

# Examples
```julia
# Rotate horizontal to 45°
R = Rotator(π/4)
j_out = R * horizontal()
```
```
"""
function Rotator(θ::Real)
    c = cos(θ)
    s = sin(θ)
    M = [c s; -s c]
    return JonesMatrix(M)
end

"""
  WavePlate(Γ::Real, θ::Real=0.0)

  General waveplate with phase retarardation Γ and fast axis at angle θ (in radians).

# Arguments
- Γ:: Real: Phase retardation in radians (2π(ne-no)d/λ)
  - `Γ = π/2`: quarter-wave plate
  - `Γ = π`: half-wave plate
- θ:: Real: Angle of the fast axis in radians
"""
function WavePlate(Γ::Real, θ::Real = 0.0)
    return Rotator(-θ) * JonesMatrix([exp(-im*Γ/2.0) 0; 0 exp(-im*Γ/2.0)]) * Rotator(θ)
end

"""

  QuarterWavePlate(θ::Real)

Quarter-wave plate (λ/4 plate) with fast axis at angle θ (in radians).

# Arguments
- θ::Real: Angle of the fast axis in radians

# Examples
```julia
# QWP with fast axis horizontal
Q = QuarterWavePlate(0)

# Convert horizontal to right circular
j_out = Q * horizontal()
```
"""
QuarterWavePlate(θ::Real) = WavePlate(π/2, θ)

""" HalfWavePlate(θ::Real)

Half-wave plate (λ/2 plate) with fast axis at angle θ (in radians).

# Arguments
- θ::Real: Angle of the fast axis in radians

# Examples
```Julia
# HWP with fast axis at 22.5° rotates horizontal to vertical
H = HalfWavePlate(π/8)
j_out = H * horizontal()
```
"""
HalfWavePlate(θ::Real) = WavePlate(π, θ)



# ============================================================================
# Analysis functions
# ============================================================================
"""
    ellipse_parameters(j::JonesVector)

Extract ellipse parameters from Jones vector.

# Returns
- `NamedTuple` with fields:
  - `a::Float64`: Semi-major axis
  - `b::Float64`: Semi-minor axis
  - `ψ::Float64`: Orientation angle (radians)
  - `χ::Float64`: Ellipticity angle (radians, -π/4 ≤ χ ≤ π/4)

# Notes
- χ > 0: right-handed ellipse
- χ < 0: left-handed ellipse
- χ = ±π/4: circular polarization
- χ = 0: linear polarization
"""
function ellipse_parameters(j::JonesVector)
    Ex, Ey = j.Ex, j.Ey

    # Normalize
    I = intensity(j)
    Ex_norm = Ex / sqrt(I)
    Ey_norm = Ey / sqrt(I)

    # Extract amplitude and phase
    ax = abs(Ex_norm)
    ay = abs(Ey_norm)
    δ = angle(Ey_norm) - angle(Ex_norm)  # Phase difference

    # Orientation angle
    ψ = 0.5 * atan(2 * ax * ay * cos(δ), ax^2 - ay^2)

    # Ellipticity angle
    χ = 0.5 * asin(clamp(2 * ax * ay * sin(δ), -1, 1))

    # Semi-axes
    a = sqrt((ax^2 + ay^2 + sqrt((ax^2 - ay^2)^2 + 4 * ax^2 * ay^2 * cos(δ)^2)) / 2)
    b = sqrt((ax^2 + ay^2 - sqrt((ax^2 - ay^2)^2 + 4 * ax^2 * ay^2 * cos(δ)^2)) / 2)

    return (a = a, b = b, ψ = ψ, χ = χ)
end


""" degree_of_polarization(j::JonesVector)

Calculate the degree of polarization (always 1.0 for pure Jones vectors).

For partially polarized light, use Stokes parameters or Mueller calculus instead. """
degree_of_polarization(j::JonesVector) = 1.0


