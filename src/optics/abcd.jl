using LinearAlgebra

export AbstractOpticalElement, FreeSpace, ThinLens, ThickLens, Interface, PlanoConvexLens
export transfer_matrix, effective_focal_length, back_focal_length

abstract type AbstractOpticalElement end

"""
    FreeSpace(L)
"""
struct FreeSpace <: AbstractOpticalElement
    L::Real
end

"""
    ThinLens(f)
Thin lens with focal length f
"""
struct ThinLens <: AbstractOpticalElement
    f::Real
end

"""
    Interface(n1, n2, R=Inf)

Interfarce between the medium with Refractive n1 and n2, R is radius of curvature.
"""
struct Interface <: AbstractOpticalElement
    n1::Real
    n2::Real
    R::Real
end

"""
    ThickLens(n_lens, n_env, d, R1, R2)

Thickness: d, Refractive Index of lens: n_lens, Refractive Index of environment: n_env, curvature radius: R1, R2
"""
struct ThickLens <: AbstractOpticalElement
    n_lens::Real
    n_env::Real
    d::Real
    R1::Real
    R2::Real
end
# --- (Transfer Matrix) ---

function transfer_matrix(fs::FreeSpace)
    return [
        1.0 fs.L;
        0.0 1.0
    ]
end

function transfer_matrix(tl::ThinLens)
    return [
        1.0 0.0;
        -1.0/tl.f 1.0
    ]
end

function transfer_matrix(it::Interface)
    P = (it.n2 - it.n1) / it.R
    return [
        1.0 0.0;
        (it.n2 - it.n1)/(it.n2*it.R) it.n1/it.n2
    ]
end

"""
    transfer_matrix(tl::ThickLens)

Calculate the ABCD transfer matrix for a thick lens.
- `tl`: ThickLens object

The transfer matrix is computed as the product of:
1. Interface from environment to lens (radius R1)
2. Free space propagation through the lens (thickness d)
3. Interface from lens to environment (radius -R2)
"""
function transfer_matrix(tl::ThickLens)
    m1 = Interface(tl.n_env, tl.n_lens, tl.R1) |> transfer_matrix
    m2 = FreeSpace(tl.d) |> transfer_matrix
    m3 = Interface(tl.n_lens, tl.n_env, -tl.R2) |> transfer_matrix
    return m3 * m2 * m1
end

"""
  ThickLens(n_lens, d, R1, R2; n_env=1.0)

Convenience constructor for ThickLens.
- `n_lens`: Refractive index of the lens
- `d`: Thickness of the lens
- `R1`: Radius of curvature for the first surface
- `R2`: Radius of curvature for the second surface
- `n_env`: (optional) Refractive index of the environment (default: 1.0)
"""
function ThickLens(n_lens, d, R1, R2; n_env = 1.0)
    return ThickLens(n_lens, n_env, d, R1, R2)
end

"""
    PlanoConvexLens(n_lens, d, R; n_env=1.0)

Convenience constructor for a plano-convex lens (one flat side, one curved side).

- `n_lens`: Refractive index of the lens
- `d`: Thickness of the lens
- `R`: Radius of curvature for the curved side
- `n_env`: (optional) Refractive index of the environment (default: 1.0)
"""
function PlanoConvexLens(n_lens, d, R; n_env = 1.0)
    return ThickLens(n_lens, n_env, d, Inf, R)
end

"""
    effective_focal_length(M)
"""
function effective_focal_length(M::AbstractMatrix)
    # 1/f = -C
    return -1.0 / M[2, 1]
end

"""
    back_focal_length(M)
"""
function back_focal_length(M::AbstractMatrix)
    # BFL = (1 - A) / -C (or A/(-C) depending on convention)
    return M[1, 1] / -M[2, 1]
end

