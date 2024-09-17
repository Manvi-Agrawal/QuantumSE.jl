using QuantumSE
using Z3
using SparseArrays
using LinearAlgebra: nullspace

include("QEC_Defaults.jl")
using .QEC_Defaults

include("QEC_Pipeline.jl")
using .QEC_Pipeline

include("QEC_Helper.jl")
using .QEC_Helper

ctx = Z3.Context()

struct PrimeG <: AbstractGroup
    p::Integer
    v::Integer

    PrimeG(p, v) = new(p, mod(v, p))
end

Base.:*(x::PrimeG, y::PrimeG) = PrimeG(x.p, x.v + y.v)
Base.inv(x::PrimeG) = PrimeG(x.p, x.p - x.v)

Hamming743 = GF2.([1 1 0 1 1 0 0;1 0 1 1 0 1 0;0 1 1 1 0 0 1])
Hamming733 = GF2.([1 0 0 0 1 1 0;0 1 0 0 1 0 1;0 0 1 0 0 1 1;0 0 0 1 1 1 1])

function TannerCode(G::Vector{<:AbstractGroup}, A::Vector{<:AbstractGroup}, B::Vector{<:AbstractGroup}, HA, HAt, HB, HBt)
    ng = length(G)
    na = length(A)
    nb = length(B)

    #  Cayley complex 
    F = [((g,0,0), (a*g,0,1), (a*g*b,1,1), (g*b,1,0)) for g in G for a in A for b in B]

    nf = length(F)

    # Tanner code
    F_idxs = Dict([(F[j], j) for j in 1:nf])

    IaX = Vector{Int}(undef, na*nb*ng*2)
    IbX = Vector{Int}(undef, na*nb*ng*2)
    IgX = Vector{Int}(undef, na*nb*ng*2)
    IfX = Vector{Int}(undef, na*nb*ng*2)
    IaZ = Vector{Int}(undef, na*nb*ng*2)
    IbZ = Vector{Int}(undef, na*nb*ng*2)
    IgZ = Vector{Int}(undef, na*nb*ng*2)
    IfZ = Vector{Int}(undef, na*nb*ng*2)
    @inbounds @simd for jg in 1:ng
        g = G[jg]
        for jb in 1:nb
            b = B[jb]
            for ja in 1:na
                a = A[ja]
                j = (jg-1)*na*nb*2+(jb-1)*na*2+(ja-1)*2+1
                IaX[j], IaX[j+1], IaZ[j], IaZ[j+1] = ja, ja, ja, ja
                IbX[j], IbX[j+1], IbZ[j], IbZ[j+1] = jb, jb, jb, jb
                IgX[j], IgX[j+1], IgZ[j], IgZ[j+1] = jg, jg+ng, jg, jg+ng
                IfX[j], IfX[j+1], IfZ[j], IfZ[j+1] = F_idxs[((g,0,0), (a*g,0,1), (a*g*b,1,1), (g*b,1,0))], F_idxs[((inv(a)*g*inv(b),0,0), (g*inv(b),0,1), (g,1,1), (inv(a)*g,1,0))], F_idxs[((inv(a)*g,0,0), (g,0,1), (g*b,1,1), (inv(a)*g*b,1,0))], F_idxs[((g*inv(b),0,0), (a*g*inv(b),0,1), (a*g,1,1), (g,1,0))]
            end
        end
    end

    IX = @. (IgX-1)*nf + IfX
    JX = @. (IaX-1)*nb + IbX
    VX = ones(GF2, na*nb*ng*2)

    IZ = @. (IgZ-1)*nf + IfZ
    JZ = @. (IaZ-1)*nb + IbZ
    VZ = ones(GF2, na*nb*ng*2)

    HXt = sparse(IX, JX, VX, nf*2*ng, na*nb)
    HZt = sparse(IZ, JZ, VZ, nf*2*ng, na*nb)

    HXt = reshape(HXt * kron(sparse(transpose(HA)), sparse(transpose(HB))), nf, :)
    HZt = reshape(HZt * kron(sparse(transpose(HAt)), sparse(transpose(HBt))), nf, :)

    HXt, HZt
end


function get_adj(HXt, HZt, nx, nz)
    X_idxs = [findall(!iszero, HXt[:,j]) for j in 1:nx]
    Z_idxs = [findall(!iszero, HZt[:,j]) for j in 1:nz]
    _xadj(i) = X_idxs[i]
    _zadj(i) = Z_idxs[i]

    return (_xadj, _zadj)
end

function get_tanner_code(m, k)
    #=
    p = 13
    q = 5

    G = pgl2(q)
    S = jacobi4squares(p, q)
    A = S
    B = S

    HA = [Hamming743 Hamming743]
    HAt = [Hamming733 Hamming733]
    HB = [Hamming733 Hamming733]
    HBt = [Hamming743 Hamming743]
    =#

    #G = [PrimeG(7^m*2^k, j-1) for j in 1:7^m*2^k]
    #A = [PrimeG(7^m*2^k, (j-1)*7^(m-1)*2^k) for j in 1:7]
    #B = [PrimeG(7^m*2^k, (j-1)*7^(m-1)*2^k) for j in 1:7]

    aa = k
    G = [PrimeG(7^m*aa, j-1) for j in 1:7^m*aa]
    A = [PrimeG(7^m*aa, (j-1)*7^(m-1)*aa) for j in 1:7]
    B = [PrimeG(7^m*aa, (j-1)*7^(m-1)*aa) for j in 1:7]
    
    #HA = Hamming313
    #HAt = Hamming322
    #HB = Hamming322
    #HBt = Hamming313

    HA = Hamming743
    HAt = Hamming733
    HB = Hamming733
    HBt = Hamming743

    return TannerCode(G, A, B, HA, HAt, HB, HBt)

end

function get_stabilizer(HXt, HZt)
    stabilizer, dx, dz = stabilizer_from_css_code(Matrix{GF2}(transpose(HXt)), Matrix{GF2}(transpose(HZt)), ctx)
    return stabilizer
end

function get_phases(HXt, HZt)
    phases, dx, dz = phases_from_css_code(Matrix{GF2}(transpose(HXt)), Matrix{GF2}(transpose(HZt)), ctx)
    return phases
end

function get_tanner_decoder_config(k::Integer)
    (HXt, HZt) = get_tanner_code(1, k)
    
    n, nx = size(HXt)
    nz = size(HZt,2)
    @show n, k, nx, nz

    (_xadj, _zadj) = get_adj(HXt, HZt, nx, nz)

    tanner_decoder_config = QEC_Pipeline.QecPipelineConfig(
        d=6, # min(dx, dz)??
        num_qubits = 343*k, # TODO: Fix this nq, currently hardcode for this specific Hamming code
        _xadj=_xadj,
        _zadj=_zadj,
        nx=nx,
        nz=nz,
        stabilizer=get_stabilizer(HXt, HZt),
        phases= get_phases(HXt, HZt),
        ctx=ctx)

    return tanner_decoder_config
end

QEC_Pipeline.qec_runner("tanner_code.csv", get_tanner_decoder_config, 1:4)
