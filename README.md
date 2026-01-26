# LaserToolBox

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/julia-v1.12-blue.svg)](https://julialang.org/)

LaserToolBox is a comprehensive suite of tools designed for laser spectroscopy system design and optimization. It provides numerical models for dispersion, refractive indices, optical propagation, and polarization analysis.

## Features

### 📐 Dispersion Analysis

- **Material Refractive Indices**: Accurate dispersion models for various optical materials
  - **Glass materials**: BK7, Fused Silica (UVFS), CaF₂, SF10, SF11
  - **Biaxial crystals**: Calcite, Quartz, MgF₂
  - **Nonlinear crystals**: α-BBO, β-BBO
  - **Air**: Temperature and pressure-dependent refractive index (Ciddor model)
- **Dispersion Models**:
  - Sellmeier equation with symbolic and numeric implementations
  - Group Velocity Dispersion (GVD)
  - Third Order Dispersion (TOD)
  - Higher-order dispersion coefficients (β_n)

- **Dispersive Optics**:
  - Prism pair calculations
  - Group Delay Dispersion (GDD)
  - Brewster angle calculations
  - Ideal apex angle optimization

### 🔬 Optical Systems (ABCD Matrix)

- Matrix-based ray transfer analysis
- Supported optical elements:
  - Free space propagation
  - Thin lenses
  - Thick lenses
  - Interfaces (refraction/reflection)
  - Plano-convex lenses
- Effective and back focal length calculations

### 🌊 Polarization Analysis (Jones Calculus)

- **Jones Vectors**: Complete polarization state representation
  - Linear: horizontal, vertical, diagonal, antidiagonal
  - Circular: left and right circular polarization
- **Jones Matrices**: Optical element modeling
  - Linear polarizers
  - Wave plates (quarter-wave, half-wave, arbitrary retardance)
  - Rotators
- **Polarization Analysis**:
  - Intensity calculations
  - Degree of polarization
  - Ellipse parameters

### 🧮 Numerical Tools

- Automatic differentiation for dispersion derivatives
- nth-order derivative calculations using ForwardDiff

## Installation

```julia
using Pkg
Pkg.add("LaserToolBox")
```

Or in the Julia REPL package mode (press `]`):

```julia
pkg> add LaserToolBox
```

## Quick Start

```julia
using LaserToolBox

# Calculate refractive index of BK7 at 800 nm
λ = 0.8  # wavelength in μm
n_bk7 = n.bk7(λ)

# Get the ordinary and extraordinary refractive indices of calcite
n_calcite = n.calcite(1.064)  # at 1064 nm
println("n_e = $(n_calcite.n_e), n_o = $(n_calcite.n_o)")

# Calculate GVD of fused silica
gvd_fs = gvd(n.fused_silica, 0.8)  # in fs²/mm

# ABCD matrix for a thin lens
lens = ThinLens(f = 100.0)  # 100 mm focal length
M = transfer_matrix(lens)

# Jones vector for 45° linear polarization
pol = diagonal()
```

## Usage Examples

### Dispersion Calculations

```julia
# Refractive index and its derivatives
λ = 0.8  # 800 nm
n_val = n.bk7(λ)                    # refractive index
dn_dλ = n.bk7(λ; derivative = 1)    # first derivative
d2n_dλ2 = n.bk7(λ; derivative = 2)  # second derivative

# Group velocity dispersion
gvd_bk7 = gvd(n.bk7, λ)  # fs²/mm
tod_bk7 = tod(n.bk7, λ)  # fs³/mm

# Prism pair for dispersion compensation
prism = PrismPair(
    material = n.sf11,
    apex_angle_deg = 68.0,
    separation_mm = 100.0
)
gdd_val = gdd(prism, λ)
```

### Optical System Design

```julia
# Design a simple lens system
lens1 = ThinLens(f = 50.0)
space = FreeSpace(L = 100.0)
lens2 = ThinLens(f = 75.0)

# Calculate total transfer matrix
M_total = transfer_matrix(lens2) * transfer_matrix(space) * transfer_matrix(lens1)

# Calculate focal length
f_eff = effective_focal_length(M_total)
f_back = back_focal_length(M_total)

# Plano-convex lens example (Edmund Optics #11-720)
pcx = PlanoConvexLens(
    n_lens = n.fused_silica(0.5876),
    d = 1.6,        # thickness in mm
    R = 114.75      # radius of curvature in mm
)
```

### Polarization Manipulation

```julia
# Create a polarization state
input_pol = horizontal()

# Quarter-wave plate at 45°
qwp = QuarterWavePlate(θ = 45.0)
output_pol = qwp * input_pol

# Analyze the output
I = intensity(output_pol)
dop = degree_of_polarization(output_pol)
ellipse = ellipse_parameters(output_pol)
```

## API Overview

### Refractive Index Module (`n` or `RefractiveIndex`)

All material functions follow the same interface:

```julia
n.material(λ; derivative=0)
```

Available materials:

- `n.air(λ)` - Air (Ciddor model)
- `n.bk7(λ)` - Schott BK7
- `n.fused_silica(λ)` - Fused silica (UVFS)
- `n.caf2(λ)` - Calcium fluoride
- `n.sf10(λ)`, `n.sf11(λ)` - Schott SF glasses
- `n.calcite(λ)`, `n.quartz(λ)`, `n.mgf2(λ)` - Biaxial crystals (returns `(n_e, n_o)`)
- `n.alpha_bbo(λ)`, `n.beta_bbo(λ)` - Beta-barium borate

### Dispersion Functions

```julia
gvd(material_func, λ)    # Group velocity dispersion [fs²/mm]
tod(material_func, λ)    # Third order dispersion [fs³/mm]
beta_n(material_func, λ; order=2)  # nth-order dispersion coefficient
```

### ABCD Matrix Optics

```julia
Medium(n1, n2)                 # Refractive interface
FreeSpace(L)                   # Free space propagation
ThinLens(f)                    # Thin lens
ThickLens(n, R1, R2, d)        # Thick lens
Interface(n1, n2, R)           # Curved interface
PlanoConvexLens(n_lens, d, R)  # Plano-convex lens

transfer_matrix(element)       # Get ABCD matrix
effective_focal_length(M)      # Calculate effective focal length
back_focal_length(M)           # Calculate back focal length
```

### Jones Calculus

**Polarization states:**

```julia
horizontal(), vertical(), diagonal(), antidiagonal()
right_circular(), left_circular()
```

**Optical elements:**

```julia
LinearPolarizer(θ)              # Linear polarizer at angle θ
QuarterWavePlate(θ)             # QWP with fast axis at θ
HalfWavePlate(θ)                # HWP with fast axis at θ
WavePlate(δ, θ)                 # Arbitrary wave plate
Rotator(θ)                      # Coordinate rotation
```

**Analysis:**

```julia
intensity(jones_vector)
degree_of_polarization(jones_vector)
ellipse_parameters(jones_vector)
```

## Documentation

For detailed documentation, see the [docs](docs/) directory or visit the documentation website (coming soon).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Ryuichi Arafune**

## Citation

If you use LaserToolBox in your research, please cite:

```bibtex
@software{lasertoolbox,
  author = {Arafune, Ryuichi},
  title = {LaserToolBox: A Julia package for laser spectroscopy system design},
  year = {2026},
  url = {https://github.com/arafune/LaserToolBox}
}
```
