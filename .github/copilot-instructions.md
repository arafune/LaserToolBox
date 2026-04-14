# LaserToolBox Copilot Instructions

## Build and test commands

Use Julia 1.12 with the project environment at the repository root.

```bash
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
```

Run the full test suite with:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Run a single test file by loading the package in the project environment and including that file directly. For example:

```bash
julia --project -e 'using Test, LaserToolBox; include("test/optics/test_abcd.jl")'
```

The repository also has a local coverage entrypoint that wraps `test/runtests.jl` and writes `lcov.info` plus `coverage_html/`:

```bash
julia --project --code-coverage test/runtests_local.jl
```

## High-level architecture

`src/LaserToolBox.jl` is the package entrypoint. It assembles the public API by including four main areas:

- `numerics/`: derivative helpers used by dispersion code
- `dispersion/`: symbolic models, material constants, numeric refractive-index functions, dispersion orders, and prism-pair calculations
- `optics/`: ABCD-matrix elements and focal-length helpers
- `polarization/`: Jones-vector and Jones-matrix types plus optical elements and analysis helpers

The dispersion stack is layered rather than flat:

- `src/dispersion/models/` defines reusable symbolic and numeric model formulas such as Sellmeier and Ciddor
- `src/dispersion/materials/constants/` stores coefficient tables and source metadata
- `src/dispersion/materials/numeric/` turns those constants into user-facing refractive-index functions in the `RefractiveIndex` module
- `src/dispersion/orders.jl` computes `beta_n`, `gvd`, and `tod` generically from any refractive-index function that supports the repository’s derivative interface
- `src/dispersion/dispersive_optics/prism_pair.jl` builds on that same interface to compute prism-pair geometry and GDD

At the top level, users are expected to access refractive-index functions through `LaserToolBox.n`, which is an alias for the exported `RefractiveIndex` module.

Tests mirror the package structure. `test/runtests.jl` is the aggregator, and subsystem tests live under matching folders such as `test/dispersion/...`, `test/optics/`, and `test/polarization/`.

## Key conventions

Dispersion code assumes wavelengths are in micrometers. Material functions, `beta_n`, `gvd`, and `tod` all use `λ` in `μm`. `PrismPair` is the exception at the API boundary: its constructor accepts angles in degrees and lengths in millimeters, converts them internally to radians and `μm`, and also converts wavelengths above `100` as if they were given in `nm`.

Material functions are designed to be composable. Numeric refractive-index functions must follow the shape `material(λ; derivative=0)` so they can be passed directly into `beta_n`, `gvd`, `tod`, and prism-pair routines. When values are outside a model’s valid wavelength range, the code raises `ArgumentError` instead of clamping or returning sentinel values.

Aggregate crystal helpers return named tuples rather than custom structs. Functions such as `calcite(λ)`, `mgf2(λ)`, `quartz(λ)`, `alpha_bbo(λ)`, and `beta_bbo(λ)` return `(n_e = ..., n_o = ...)`, while `_e` and `_o` helper functions expose each branch separately.

Optics and polarization use different angle conventions. Jones-calculus constructors such as `LinearPolarizer`, `Rotator`, and wave plates take angles in radians. Prism-pair helpers with `_deg` in the name return degrees, and the `PrismPair` constructor accepts degree inputs before internal conversion.

Formatting uses JuliaFormatter’s `blue` style via `.JuliaFormatter.toml`.
