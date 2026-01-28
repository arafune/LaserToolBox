"""
Import the refractive index function for SF11 glass from the Materials.RefractiveIndex module.
This function is used as the default material property for the PrismPair struct.
"""
# dispersion/dispersion_optics/prism_pair.jl
#
# PrismPair to calculate GDD
using ..Materials.RefractiveIndex: sf11
using ...Dispersion: beta_n, gvd

"""

  A struct representing a pair of prisms used for dispersion compensation.

# Fields
- `incident_angle::Float64`  
  Incident angle of the beam (radians; input in degrees, converted internally).
- `material::F`  
  Function to get refractive index: n(λ) (λ in μm). Must accept a `derivative` keyword argument.
- `apex_angle::Float64`  
  Apex angle of the prism (radians; input in degrees, converted internally).
- `separation::AbstractVector{Float64}`
  Separation between the prisms (μm; input in mm, converted internally).
- `insertion::Tuple{Float64,Float64}`  
  Insertion depth of the beam into each prism (μm; input in mm, converted internally).
- `wavelength::AbstractVector{Float64}`  
Central wavelength(s) (μm; converted internally. While nm can be accepted, but strongly recommend to use explicit μm) .

# Usage
Create a `PrismPair` using keyword arguments. All units are automatically converted as described above.

# Example
```julia
pp = PrismPair(
    incident_angle = 60.0,      # degrees
    material = n.sf11,            # function n(λ)
    apex_angle = 59.0,          # degrees
    separation = 100.0,         # mm
    insertion = (3.0, 5.0),    # mm
    wavelength = .800          # nm or μm
)
```

# Note

The unit of each field is as follows:
- angle: radians
- length : µm
- time : fs

The unit of args are converted as follows:

* wavelength: µm  → converted to µm, internally 
* separation: mm  →  converted to µm, internally.
* insertion: mm → converted to µm, internally.
* apex_angle : degrees converted to radians, internally.
* incident_angle : degrees converted to radians, internally.
* material : function that takes wavelength (in µm) and returns refractive index.:

  The argument material must be a function that takes the wavelength λ as its first
  argument and supports a derivative keyword argument, e.g., material(λ; derivative=1).

"""
struct PrismPair{F}
    incident_angle::Float64 # Incident angle of the beam in radians
    material::F  # Function to get refractive index: n(λ)
    apex_angle::Float64 # Apex angle of the prism in radians
    separation::AbstractVector{Float64} # Separation between the prisms
    insertion::Tuple{Float64,Float64}  # Insertion depth of the beam into the prism
    wavelength::AbstractVector{Float64} # Central wavelength in µm
end

nm_to_μm(λ::Real) = λ > 100 ? λ * 1e-3 : λ
nm_to_μm(λ::Union{AbstractVector,Tuple}) = nm_to_μm.(λ)

# document is given in the struc
function PrismPair(;
    incident_angle = 60.0,
    material = sf11,
    apex_angle = 59.0,
    separation = 100.0,
    insertion = (3, 5),
    wavelength = 800.0,
)

    wavelength =
        wavelength isa AbstractVector ? Vector{Float64}(wavelength) : [float(wavelength)]
    wavelength = nm_to_μm(wavelength)
    separation =
        separation isa AbstractVector ? Vector{Float64}(separation * 1e3) :
        [float(separation*1e3)]  # convert mm to µm
    insertion = (insertion[1] * 1e3, insertion[2] * 1e3)  # convert mm to µm
    incident_angle = deg2rad(incident_angle)  # convert degrees to radians
    apex_angle = deg2rad(apex_angle)  # convert degrees to radians

    return PrismPair(
        incident_angle,
        material,
        apex_angle,
        separation,
        insertion,
        wavelength,
    )
end



"""
theta_1(prism_pair::PrismPair)

The angle of refraction inside the first prism.

Returns a scalar or vector depending on wavelength.
"""
function theta_1(prism_pair::PrismPair)
    n = prism_pair.material.(prism_pair.wavelength)
    θ0 = prism_pair.incident_angle
    return asin.(sin(θ0) ./ n)
end

"""
theta_2(prism_pair::PrismPair)

The angle of exit outside the first prism.

Returns a scalar or vector depending on wavelength.
"""
theta_2(prism_pair::PrismPair) = begin
    θ1 = theta_1(prism_pair)
    return @. prism_pair.apex_angle - θ1
end

"""
lg(prism_pair::PrismPair)

Calculate the geometric path length through the prism pair.

Returns a scalar or vector depending on wavelength.
"""
function lg(prism_pair::PrismPair)
    θ1 = theta_1(prism_pair)
    θ2 = theta_2(prism_pair)

    l2 = prism_pair.insertion[2]
    l1 = @. prism_pair.insertion[1] * cos(θ1) / cos(θ2)
    g = @. (l1 + l2) * sin(prism_pair.apex_angle)
    return @. g / cos(θ1)
end

function dn_dΩ(prism_pair::PrismPair)
    c = 0.299792458  # µm/fs
    dn_dλ = prism_pair.material.(prism_pair.wavelength; derivative = 1)
    return @. -(prism_pair.wavelength^2) / (2pi * c) * dn_dλ
end

function dθ1_dΩ(prism_pair::PrismPair)
    n = prism_pair.material.(prism_pair.wavelength)
    θ0 = prism_pair.incident_angle
    dn_dΩ_val = dn_dΩ(prism_pair)
    return @. -(n^2 - sin(θ0)^2)^(-1 / 2) * sin(θ0) * dn_dΩ_val / n
end

dθ3_dΩ(prism_pair::PrismPair) = begin
    n = prism_pair.material.(prism_pair.wavelength)
    θ2 = theta_2(prism_pair)
    dθ1_dΩ_val = dθ1_dΩ(prism_pair)
    dn_dΩ_val = dn_dΩ(prism_pair)
    return @. (1.0 - n^2 * sin(θ2)^2)^(-1 / 2) *
              (-n * cos(θ2) * dθ1_dΩ_val + sin(θ2) * dn_dΩ_val)
end

raw"""
  `gdd_positive`(prism_pair::PrismPair) -> Vector{Float64}

Returns a Vector{Float64} even though wavelength was a single value.

This positive contribution to GDD arises from the material dispersion of the prism material over distance Lg.
"""
gdd_positive(prism_pair::PrismPair) = begin
    lg_val = lg(prism_pair)
    return @. gvd(prism_pair.material, prism_pair.wavelength; unit = :μm) * lg_val
end

raw"""
  `gdd_negative`(prism_pair::PrismPair) -> Vector{Float64}

Returns a Vector{Float64} even though wavelength was a single value.

This negative contribution to GDD arises from the angular dispersion $\frac{d\theta_3}{d\Omega},
and angular dispersion of $\frac{d\theta_3}{d\Omega}$ over distance $L_g$ in the glass of refractive index n.
"""
gdd_negative(prism_pair::PrismPair) = begin
    dn_dλ = prism_pair.material.(prism_pair.wavelength; derivative = 1)
    n = prism_pair.material.(prism_pair.wavelength)
    dθ3_dΩ_val = dθ3_dΩ(prism_pair)
    dθ1_dΩ_val = dθ1_dΩ(prism_pair)
    lg_val = lg(prism_pair)
    gap_component =
        @. -2pi / prism_pair.wavelength * prism_pair.separation * dθ3_dΩ_val^2
    prism_component = @. -n * 2pi / prism_pair.wavelength * lg_val * (dθ1_dΩ_val)^2
    return @. (gap_component + prism_component)
end


"""
  `gdd`(prism_pair::PrismPair) -> Float64 or Vector{Float64}

GDD of Prism pair

Returns a scalar Float64 if wavelength was a single value, otherwise returns a Vector{Float64}.

# Note

The total GDD is the sum of the positive GDD from material dispersion and the negative GDD from geometric dispersion.
This function returns the GDD for a single prism pair by summing both contributions.
In the actual setup, the beam passes through two prism pairs, so the total GDD is **twice** the value returned by this function.

# Reference

 Ultrashort Laser Pulse Phenomena (2006).

 W. R. Jean-Claude Diels
"""
gdd(prism_pair::PrismPair) = begin
    result = gdd_positive(prism_pair) .+ gdd_negative(prism_pair)
    # Return scalar for single-wavelength input (backward compatibility)
    return length(result) == 1 ? result[1] : result
end

"""
  `brewster_angle_deg`(material::Function, λ)

Return the Brewster angle for a given material at λ.

# Arguments
- `material::Function, λ`: A function describing the refractive index as a function of wavelength `λ`. Returns the Brewster angle for each `λ`.
- `n::Real`: A real refractive index. Returns the Brewster angle for this value.
- `n::AbstractArray`: An array of refractive indices. Returns an array of Brewster angles for each value.

# Example
```julia
brewster_angle_deg(n)             # n is a real number
brewster_angle_deg([1.5, 2.0])    # n is an array
brewster_angle_deg(material, λ)   # material is a function of λ
```
"""
function brewster_angle_deg(material::Function, λ)
    return @. rad2deg(atan(material(λ)))
end

"""
  `ideal_apex_deg`(material::Function, λ::Float64) -> Float64

Return the ideal apex angle for the prism at λ.

# Arguments
- `material::Function, λ`: A function describing the refractive index as a function of wavelength `λ`. Returns the Brewster angle for each `λ`.
- `n::Real`: A real refractive index. Returns the Brewster angle for this value.
- `n::AbstractArray`: An array of refractive indices. Returns an array of Brewster angles for each value.
"""
function ideal_apex_deg(material::Function, λ)
    return @. rad2deg(2 * asin(sin(deg2rad(brewster_angle_deg(material, λ)))/material(λ)))
end
