using Test
using LaserToolBox
using LinearAlgebra

@testset "JonesVector" begin
    @testset "Construction" begin
        # From complex numbers
        j1 = JonesVector(1.0 + 0.0im, 0.0 + 0.0im)
        @test j1.Ex == 1.0 + 0.0im
        @test j1.Ey == 0.0 + 0.0im
        
        # From real numbers
        j2 = JonesVector(1.0, 0.0)
        @test j2.Ex == 1.0 + 0.0im
        @test j2.Ey == 0.0 + 0.0im
    end
    
    @testset "Vector conversion" begin
        j = JonesVector(1.0, 2.0im)
        v = Vector(j)
        @test v == [1.0 + 0.0im, 0.0 + 2.0im]
        @test Array(j) == v
    end    
    @testset "Arithmetic operations" begin
        j1 = JonesVector(1.0, 2.0)
        j2 = JonesVector(3.0, 4.0)
        
        # Addition
        j_sum = j1 + j2
        @test j_sum.Ex == 4.0
        @test j_sum.Ey == 6.0
        
        # Subtraction
        j_diff = j2 - j1
        @test j_diff.Ex == 2.0
        @test j_diff.Ey == 2.0
        
        # Scalar multiplication
        j_mul = 2.0 * j1
        @test j_mul.Ex == 2.0
        @test j_mul.Ey == 4.0
        
        j_mul2 = j1 * 2.0
        @test j_mul2.Ex == 2.0
        @test j_mul2.Ey == 4.0
        
        # Division
        j_div = j1 / 2.0
        @test j_div.Ex == 0.5
        @test j_div.Ey == 1.0
    end
    
    @testset "Norm and normalization" begin
        j = JonesVector(3.0, 4.0)
        @test norm(j) ≈ 5.0
        
        j_norm = normalize(j)
        @test norm(j_norm) ≈ 1.0
        @test j_norm.Ex ≈ 0.6
        @test j_norm.Ey ≈ 0.8
    end
    
    @testset "Intensity" begin
        j = JonesVector(3.0, 4.0)
        @test intensity(j) ≈ 25.0
        
        j2 = JonesVector(1.0, 1.0im)
        @test intensity(j2) ≈ 2.0
    end
end

@testset "JonesMatrix" begin
    @testset "Construction" begin
        # From complex matrix
        M_complex = ComplexF64[1 0; 0 1]
        J1 = JonesMatrix(M_complex)
        @test J1.M == M_complex
        
        # From real matrix
        M_real = [1.0 0.0; 0.0 1.0]
        J2 = JonesMatrix(M_real)
        @test J2.M == ComplexF64.(M_real)
    end
    
    @testset "Matrix conversion" begin
        M = ComplexF64[1 2; 3 4]
        J = JonesMatrix(M)
        @test Matrix(J) == M
        @test Array(J) == M
    end
    
    @testset "Matrix-vector multiplication" begin
        # Identity matrix
        J_id = JonesMatrix([1.0 0.0; 0.0 1.0])
        j = JonesVector(1.0, 2.0)
        j_result = J_id * j
        @test j_result.Ex ≈ 1.0
        @test j_result.Ey ≈ 2.0
        
        # Horizontal polarizer
        J_h = JonesMatrix([1.0 0.0; 0.0 0.0])
        j_vert = JonesVector(0.0, 1.0)
        j_blocked = J_h * j_vert
        @test j_blocked.Ex ≈ 0.0
        @test j_blocked.Ey ≈ 0.0
    end
    
    @testset "Matrix-matrix multiplication" begin
        J1 = JonesMatrix([1.0 2.0; 3.0 4.0])
        J2 = JonesMatrix([5.0 6.0; 7.0 8.0])
        J_prod = J1 * J2
        expected = [1.0 2.0; 3.0 4.0] * [5.0 6.0; 7.0 8.0]
        @test Matrix(J_prod) ≈ expected
    end
    
    @testset "Inverse and determinant" begin
        J = JonesMatrix([2.0 1.0; 1.0 2.0])
        J_inv = inv(J)
        J_prod = J * J_inv
        @test Matrix(J_prod) ≈ [1.0 0.0; 0.0 1.0] atol=1e-10
        
        @test det(J) ≈ 3.0
    end
end

@testset "Common polarization states" begin
    @testset "Linear polarizations" begin
        h = horizontal()
        @test h.Ex ≈ 1.0
        @test h.Ey ≈ 0.0
        @test norm(h) ≈ 1.0
        
        v = vertical()
        @test v.Ex ≈ 0.0
        @test v.Ey ≈ 1.0
        @test norm(v) ≈ 1.0
        
        d = diagonal()
        @test d.Ex ≈ 1/sqrt(2)
        @test d.Ey ≈ 1/sqrt(2)
        @test norm(d) ≈ 1.0
        
        a = antidiagonal()
        @test a.Ex ≈ 1/sqrt(2)
        @test a.Ey ≈ -1/sqrt(2)
        @test norm(a) ≈ 1.0
    end
    
    @testset "Circular polarizations" begin
        r = right_circular()
        @test abs(r.Ex) ≈ 1/sqrt(2)
        @test abs(r.Ey) ≈ 1/sqrt(2)
        @test angle(r.Ey) - angle(r.Ex) ≈ π/2
        @test norm(r) ≈ 1.0
        
        l = left_circular()
        @test abs(l.Ex) ≈ 1/sqrt(2)
        @test abs(l.Ey) ≈ 1/sqrt(2)
        @test angle(l.Ey) - angle(l.Ex) ≈ -π/2
        @test norm(l) ≈ 1.0
    end
end

@testset "Optical elements" begin
    @testset "LinearPolarizer" begin
        # Horizontal polarizer
        P_h = LinearPolarizer(0.0)
        h = horizontal()
        v = vertical()
        
        j_h = P_h * h
        @test intensity(j_h) ≈ 1.0
        
        j_v = P_h * v
        @test intensity(j_v) ≈ 0.0 atol=1e-10
        
        # 45° polarizer
        P_45 = LinearPolarizer(π/4)
        d = diagonal()
        j_45 = P_45 * d
        @test intensity(j_45) ≈ 1.0 atol=1e-10
    end
    
    @testset "QuarterWavePlate" begin
        # QWP with fast axis horizontal converts diagonal to circular
        Q = QuarterWavePlate(0.0)
        h = horizontal()
        
        # Horizontal through QWP stays horizontal
        j_out = Q * h
        @test abs(j_out.Ex) ≈ 1.0
        @test abs(j_out.Ey) ≈ 0.0 atol=1e-10
    end
    
    @testset "HalfWavePlate" begin
        # HWP with fast axis at 0° rotates polarization
        H = HalfWavePlate(0.0)
        h = horizontal()
        j_out = H * h
        @test abs(j_out.Ex) ≈ 1.0
        
        v = vertical()
        j_v = H * v
        @test abs(j_v.Ey) ≈ 1.0
    end
    
    @testset "Rotator" begin
        # Rotate horizontal to vertical
        R = Rotator(π/2)
        h = horizontal()
        j_rot = R * h
        @test abs(j_rot.Ex) ≈ 0.0 atol=1e-10
        @test abs(j_rot.Ey) ≈ 1.0
        
        # Rotate to 45°
        R_45 = Rotator(π/4)
        j_45 = R_45 * h
        @test abs(j_45.Ex) ≈ 1/sqrt(2) atol=1e-10
        @test abs(j_45.Ey) ≈ 1/sqrt(2) atol=1e-10
    end
end

@testset "Analysis functions" begin
    @testset "ellipse_parameters" begin
        # Horizontal linear polarization
        h = horizontal()
        params_h = ellipse_parameters(h)
        @test params_h.χ ≈ 0.0 atol=1e-10
        
        # Vertical linear polarization
        v = vertical()
        params_v = ellipse_parameters(v)
        @test params_v.χ ≈ 0.0 atol=1e-10
        
        # Right circular polarization
        r = right_circular()
        params_r = ellipse_parameters(r)
        @test abs(params_r.χ) ≈ π/4 atol=1e-10
        @test params_r.a ≈ params_r.b atol=1e-10
        
        # Left circular polarization
        l = left_circular()
        params_l = ellipse_parameters(l)
        @test abs(params_l.χ) ≈ π/4 atol=1e-10
        @test params_l.a ≈ params_l.b atol=1e-10
    end
    
    @testset "degree_of_polarization" begin
        j = JonesVector(1.0, 1.0im)
        @test degree_of_polarization(j) == 1.0
    end
end
