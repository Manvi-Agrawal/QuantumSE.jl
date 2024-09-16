using QuantumSE
using Z3

include("QEC_Pipeline.jl")
using .QEC_Pipeline

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

@qprog repetition_decoder (nq) begin
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

    println("Type(phases): $(typeof(phases))")
    println("Length(phases): $(length(phases))")

    return phases
end

function get_stabilizer(n)
    num_qubits = n

    stabilizer = Matrix{Bool}(undef, num_qubits, 2*num_qubits)

    @simd for i in 1:n-1
        stabilizer[i+1,:] = repetition_s(n, i)
    end

    stabilizer[1,:] = repetition_lx(n)

    println("Type(stabilizer): $(typeof(stabilizer))")

    return stabilizer
end

function get_repetition_decoder_config(n)
    repetition_decoder_params = (n)
    repetition_decoder_config = QEC_Pipeline.QecDecoderConfig(
        d=n,
        num_qubits=n,
        stabilizer=get_stabilizer(n),
        phases= get_phases(n),
        ctx=ctx,
        _zadj=_adj,
        decoder=repetition_decoder,
        decoder_params=repetition_decoder_params)

    # repetition_decoder_config.x_syndrome_circuit = repetition_m_xx
    repetition_decoder_config.z_syndrome_circuit = repetition_m_zz
    repetition_decoder_config.decoder_algo_xz = mwpm

    return repetition_decoder_config
end

QEC_Pipeline.qec_runner("repetition_code.csv", get_repetition_decoder_config, 3:5)

