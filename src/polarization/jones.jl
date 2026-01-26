function WavePlate(Γ::Real, θ::Real=0.0)
    # Wave plate aligned at 0° with symmetric phase factors
    M_plate = JonesMatrix([exp(-1im*Γ/2.0) 0.0; 0.0 exp(1im*Γ/2.0)])
    
    # Rotate to desired angle: R(-θ) * M_plate * R(θ)
    return Rotator(-θ) * M_plate * Rotator(θ)
end