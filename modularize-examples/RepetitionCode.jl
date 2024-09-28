using QuantumSE
using Z3

include("QEC_Pipeline.jl")
using .QEC_Pipeline

include("QecDecoder.jl")
include("QecPipelineConfig.jl")

ctx = Context()

_adj(idx) = [idx, idx+1]

function mwpm(n::Integer, s)

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    r = [_bv_const(ctx, "r_$(j)") for j in 1:n]
    for j in 1:n-1
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_adj(j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(n)), x)).(r) ) <= bv_val(ctx, (n-1)÷2, _len2(n)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

@qprog repetition_m_zz (n, idx) begin
    CNOT(idx, idx%n+1)
    res = M(idx%n+1)
    CNOT(idx, idx%n+1)
    
    res
end

function repetition_s(n, idx)
    s = zeros(Bool, 2*n)

    s[n+idx] = true
    s[n+idx+1] = true

    s
end

function repetition_lx(n)
    s = zeros(Bool, 2*n)

    s[n+1] = true

    s
end

function repetition_lz(n)
    s = zeros(Bool, 2*n)

    for j in 1:n
        s[j] = true
    end

    s
end

@qprog repetition_decoder_ckt (nq) begin
    s = [z_syndrome_circuit(nq, j) for j in 1:nq]

    r = decoder_algo_xz(nq, s)

    for j in 1:nq
        sX(j, r[j])
    end

    e = reduce(&, r[1:((nq-1)÷2)])

    sX(1, e)
end

function get_phases(nq)
    phases = Vector{Z3.ExprAllocated}(undef, nq)
    lx = _bv_const(ctx, "lx")
    lz = _bv_const(ctx, "lz")

    for i in 1:nq-1
        phases[i+1] = _bv_val(ctx, 0)
    end

    phases[1] = lx

    return phases
end

function get_stabilizer(nq::Integer)
    stabilizer = Matrix{Bool}(undef, nq, 2*nq)

    @simd for i in 1:nq-1
        stabilizer[i+1,:] = repetition_s(nq, i)
    end

    stabilizer[1,:] = repetition_lx(nq)
    return stabilizer
end

function get_repetition_decoder(nq)
    return QecDecoder(
        d=nq,
        num_qubits = nq,
        ctx=ctx,
        decoder = repetition_decoder_ckt,
        decoder_params = (nq),
        z_syndrome_circuit = repetition_m_zz,
        decoder_algo_xz = mwpm
    )

end

function get_repetition_decoder_config(nq)
    repetition_decoder_config = QecPipelineConfig(
        stabilizer=get_stabilizer(nq),
        phases= get_phases(nq),
        decoder=get_repetition_decoder(nq)
    )

    return repetition_decoder_config
end

QEC_Pipeline.qec_runner("repetition_code.csv", get_repetition_decoder_config, 3:5)
