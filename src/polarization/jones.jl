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
    χ = 0.5 * asin(2 * ax * ay * sin(δ))
    
    # Semi-axes
    a = sqrt((ax^2 + ay^2 + sqrt((ax^2 - ay^2)^2 + 4 * ax^2 * ay^2 * cos(δ)^2)) / 2)
    b = sqrt((ax^2 + ay^2 - sqrt((ax^2 - ay^2)^2 + 4 * ax^2 * ay^2 * cos(δ)^2)) / 2)
    
    return (a=a, b=b, ψ=ψ, χ=χ)
end
