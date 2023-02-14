using TPLib
using Test

∞ = Infinite()

function ≡(x,y)
    return all(map(v -> typeof(v[1]) == typeof(v[2]) && v[1] == v[2], zip(x,y)))
end

tropicalidentity(n) = [i == j ? 0 : ∞ for i in 1:n, j in 1:n]
tropicalantiidentity(n) = [i+j==n ? 0 : ∞ for i in 1:n, j in 1:n]


A11 = [-Inf 2 -1 4 -Inf 2 2 -Inf 4 2; 1 -Inf 4 5 4 -Inf 4 -Inf 5 4]
A12 = [+Inf 2 -1 4 +Inf 2 2 +Inf 4 2; 1 +Inf 4 5 4 +Inf 4 +Inf 5 4]
A13 = [∞ 2 -1 4 ∞ 2 2 ∞ 4 2; 1 ∞ 4 5 4 ∞ 4 ∞ 5 4]
S11 = [0 0 0 ∞ ∞; 0 0 ∞ ∞ 0; ∞ ∞ ∞ 0 2; 0 ∞ ∞ -2 ∞; ∞ ∞ 0 ∞ -3; 0 ∞ 3 ∞ ∞; ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ 0 ∞; ∞ 0 ∞ ∞ 0; ∞ 0 ∞ -1 ∞; ∞ 0 0 ∞ ∞]
S12 = [0 ∞ ∞ ∞ -3; 0 ∞ ∞ -4 ∞; 0 -3 ∞ ∞ ∞; ∞ 0 ∞ ∞ ∞; ∞ ∞ ∞ 0 ∞; ∞ ∞ ∞ ∞ 0; ∞ ∞ 0 ∞ -3; ∞ ∞ 0 -5 ∞; ∞ 0 3 ∞ ∞]

@testset "compute_ext_rays" begin
    # Examples on empty coeffs
    @test compute_ext_rays(Matrix{Int64}(undef,0,0),0) ≡ Matrix{MaxPlus{Int64}}(undef,0,0)
    @test compute_ext_rays(Matrix{Float64}(undef,0,0),0) ≡ Matrix{MaxPlus{Float64}}(undef,0,0)
    @test compute_ext_rays(Matrix{Int64}(undef,0,0),0,:min) ≡ Matrix{MinPlus{Int64}}(undef,0,0)
    @test compute_ext_rays(Matrix{Float64}(undef,0,0),0,:min) ≡ Matrix{MinPlus{Float64}}(undef,0,0)
    @test compute_ext_rays(Matrix{Int64}(undef,0,0),5) ≡ convert(Matrix{MaxPlus{Int64}},tropicalidentity(5))
    @test compute_ext_rays(Matrix{Float64}(undef,0,0),5) ≡ convert(Matrix{MaxPlus{Float64}},tropicalidentity(5))
    @test compute_ext_rays(Matrix{Int64}(undef,0,0),5,:min) ≡ convert(Matrix{MinPlus{Int64}},tropicalidentity(5))
    @test compute_ext_rays(Matrix{Float64}(undef,0,0),5,:min) ≡ convert(Matrix{MinPlus{Float64}},tropicalidentity(5))

    # Examples A11-13 compared to S11-12
    @test compute_ext_rays(A11,5) ≡ convert(Matrix{MaxPlus{Float64}},S11)
    @test compute_ext_rays(A13,5) ≡ convert(Matrix{MaxPlus{Int64}},S11)
    @test compute_ext_rays(convert(Matrix{MaxPlus{Int64}},A13),5) ≡ convert(Matrix{MaxPlus{Int64}},S11)
    @test compute_ext_rays(convert(Matrix{MaxPlus{Float64}},A13),5) ≡ convert(Matrix{MaxPlus{Float64}},S11)
    @test compute_ext_rays(A12,5,:min) ≡ convert(Matrix{MinPlus{Float64}},S12)
    @test compute_ext_rays(A13,5,:min) ≡ convert(Matrix{MinPlus{Int64}},S12)
    @test compute_ext_rays(convert(Matrix{MinPlus{Int64}},A13),5) ≡ convert(Matrix{MinPlus{Int64}},S12)
    @test compute_ext_rays(convert(Matrix{MinPlus{Float64}},A13),5) ≡ convert(Matrix{MinPlus{Float64}},S12)
end


A21 = [-Inf 4 2 3 4; 4 -Inf -5 4 -1; -Inf -Inf 2 3 -1]
A22 = [+Inf 4 2 3 4; 4 +Inf -5 4 -1; +Inf +Inf 2 3 -1]
A23 = [∞ 4 2 3 4; 4 ∞ -5 4 -1; ∞ ∞ 2 3 -1]
S21 = [-7 ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ -2; -5 0 -3 ∞ ∞ ∞ ∞ ∞ ∞ 0; 0 ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞; -4 ∞ ∞ ∞ 0 ∞ ∞ ∞ -4 ∞; -1 ∞ 0 ∞ ∞ ∞ ∞ ∞ -1 ∞; ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ -4; ∞ ∞ 0 ∞ ∞ -9 ∞ ∞ ∞ ∞; ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞ -9 ∞; ∞ ∞ 0 -7 ∞ ∞ ∞ ∞ ∞ -2; ∞ 0 -3 -5 ∞ ∞ ∞ ∞ ∞ 0; ∞ ∞ ∞ 0 ∞ 0 ∞ ∞ ∞ ∞; ∞ ∞ ∞ ∞ 0 -5 ∞ ∞ ∞ ∞; ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ -5 ∞; ∞ ∞ 0 ∞ ∞ ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ -1 ∞ ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ 0 ∞; ∞ ∞ ∞ ∞ 0 ∞ ∞ -3 ∞ ∞; ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ 0; ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞; ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞; ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞; 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞; ∞ 0 ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞; ∞ 0 ∞ -4 ∞ ∞ ∞ ∞ ∞ 0; ∞ ∞ 0 ∞ ∞ ∞ -2 ∞ ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ -1 ∞ ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ ∞ -1; ∞ ∞ ∞ ∞ 0 ∞ 0 ∞ ∞ ∞]
S22 = [∞ ∞ 9 ∞ ∞ 0 ∞ ∞ ∞ 7; ∞ ∞ 9 ∞ ∞ 0 7 ∞ ∞ 12; 0 ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞; ∞ ∞ ∞ ∞ 5 0 ∞ ∞ 1 ∞; ∞ ∞ 9 ∞ ∞ 0 ∞ ∞ 8 ∞; 0 ∞ ∞ ∞ ∞ ∞ ∞ 9 ∞ ∞; 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞ 0 ∞; ∞ ∞ ∞ ∞ 5 ∞ ∞ ∞ 0 ∞; ∞ ∞ 9 ∞ ∞ ∞ ∞ ∞ 0 ∞; ∞ ∞ 4 ∞ ∞ ∞ ∞ ∞ ∞ 0; 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞ 5; ∞ ∞ ∞ ∞ 4 ∞ ∞ 8 0 ∞; ∞ ∞ 0 ∞ ∞ ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ ∞ 1 ∞ ∞; ∞ ∞ ∞ ∞ 3 ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ 0 ∞; ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞ 0; ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞ 0; ∞ ∞ ∞ ∞ ∞ ∞ ∞ ∞ 0 ∞; ∞ ∞ ∞ ∞ ∞ ∞ ∞ 0 ∞ ∞; ∞ ∞ ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞; ∞ ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞ ∞; ∞ 0 ∞ ∞ ∞ ∞ 0 ∞ ∞ ∞; ∞ ∞ ∞ 1 ∞ ∞ 0 ∞ ∞ 5; ∞ 0 ∞ ∞ ∞ ∞ ∞ 2 ∞ ∞; ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ 1 ∞; ∞ 0 ∞ ∞ ∞ ∞ ∞ ∞ ∞ 0; ∞ ∞ ∞ 1 ∞ ∞ ∞ ∞ ∞ 0]

@testset "compute_ext_rays_polar" begin
    # Examples on empty coeffs
    @test compute_ext_rays_polar(Matrix{Int64}(undef,0,0),0) ≡ Matrix{MaxPlus{Int64}}(undef,0,0)
    @test compute_ext_rays_polar(Matrix{Float64}(undef,0,0),0) ≡ Matrix{MaxPlus{Float64}}(undef,0,0)
    @test compute_ext_rays_polar(Matrix{Int64}(undef,0,0),0,:min) ≡ Matrix{MinPlus{Int64}}(undef,0,0)
    @test compute_ext_rays_polar(Matrix{Float64}(undef,0,0),0,:min) ≡ Matrix{MinPlus{Float64}}(undef,0,0)
    @test compute_ext_rays_polar(Matrix{Int64}(undef,0,0),5) ≡ convert(Matrix{MaxPlus{Int64}},tropicalidentity(10))
    @test compute_ext_rays_polar(Matrix{Float64}(undef,0,0),5) ≡ convert(Matrix{MaxPlus{Float64}},tropicalidentity(10))
    @test compute_ext_rays_polar(Matrix{Int64}(undef,0,0),5,:min) ≡ convert(Matrix{MinPlus{Int64}},tropicalidentity(10)[[collect(6:10);collect(1:5)],:])
    @test compute_ext_rays_polar(Matrix{Float64}(undef,0,0),5,:min) ≡ convert(Matrix{MinPlus{Float64}},tropicalidentity(10)[[collect(6:10);collect(1:5)],:])

    # Examples A21-23 compared to S21-22
    @test compute_ext_rays_polar(A21,5) ≡ convert(Matrix{MaxPlus{Float64}},S21)
    @test compute_ext_rays_polar(A23,5) ≡ convert(Matrix{MaxPlus{Int64}},S21)
    @test compute_ext_rays_polar(convert(Matrix{MaxPlus{Int64}},A23),5) ≡ convert(Matrix{MaxPlus{Int64}},S21)
    @test compute_ext_rays_polar(convert(Matrix{MaxPlus{Float64}},A23),5) ≡ convert(Matrix{MaxPlus{Float64}},S21)
    @test compute_ext_rays_polar(A22,5,:min) ≡ convert(Matrix{MinPlus{Float64}},S22)
    @test compute_ext_rays_polar(A23,5,:min) ≡ convert(Matrix{MinPlus{Int64}},S22)
    @test compute_ext_rays_polar(convert(Matrix{MinPlus{Int64}},A23),5) ≡ convert(Matrix{MinPlus{Int64}},S22)
    @test compute_ext_rays_polar(convert(Matrix{MaxPlus{Float64}},A23),5) ≡ convert(Matrix{MaxPlus{Float64}},S21)
end

A31 = [0 4 0; 0 3 4; 0 0 2]
S31 = [0 4 0; 0 4 2; 0 4 4; 0 0 2; 0 1 2]
T31 = [[3], [2, 3], [1], [2], [1, 2]]
S32 = [0 4 0; 0 0 0; 0 2 4; 0 3 0]
T32 = [[2], [1], [3], [2, 3]]

@testset "compute_halfspaces" begin
    @test compute_halfspaces(A31,3) ≡ (convert(Matrix{MaxPlus{Int64}},S31),T31)
    @test compute_halfspaces(convert(Matrix{MaxPlus{Float64}},A31),3) ≡ (convert(Matrix{MaxPlus{Float64}},S31),T31)
    @test compute_halfspaces(A31,3,:min) ≡ (convert(Matrix{MinPlus{Int64}},S32),T32)
    @test compute_halfspaces(convert(Matrix{MinPlus{Float64}},A31),3,:min) ≡ (convert(Matrix{MinPlus{Float64}},S32),T32)
end

A41 = [0 4 0; 0 3 4; 0 0 2]
S41 = [1 1 2; 1 4 2; 1 4 4; 1 3 4; 1 0 2; 1 4 0]
T41 = [[1, 2, 3, 4], [1, 5], [2, 6]]
S42 = [1 4 0; 1 3 0; 1 3 4; 1 2 4; 1 0 0; 1 0 2]
T42 = [[1, 2], [2, 3, 4, 5, 6]]    

@testset "compute_tropical_complex" begin
    @test compute_tropical_complex(A41,3) ≡ (convert(Matrix{MaxPlus{Int64}},S41),T41)
    @test compute_tropical_complex(convert(Matrix{MaxPlus{Float64}},A41),3) ≡ (convert(Matrix{MaxPlus{Float64}},S41),T41)
    @test compute_tropical_complex(A41,3,:min) ≡ (convert(Matrix{MinPlus{Int64}},S42),T42)
    @test compute_tropical_complex(convert(Matrix{MinPlus{Float64}},A41),3,:min) ≡ (convert(Matrix{MinPlus{Float64}},S42),T42)
end

A51 = [∞ 0 ∞ 1 ∞ ∞; ∞ -4 -3 0 ∞ ∞; ∞ ∞ -1 0 -6 ∞; 0 ∞ ∞ ∞ ∞ -4; 0 -8 ∞ ∞ ∞ -3]

@testset "compute_tangent_hypergraph" begin
end
