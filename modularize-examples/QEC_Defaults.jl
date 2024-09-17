module QEC_Defaults


using QuantumSE
using Z3


function decoder_algo_xz(ctx, d, s, s_type, nq, adj)

    ## pre-condition
    ϕ₁ = bool_val(ctx, true)

    ## post-condition
    ϕ₂ = bool_val(ctx, true)
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:nq]
    for j in 1:length(s)
        ϕ₂ = ϕ₂ & (s[j] ⊻ reduce(⊻,  r[adj(j)]) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(nq)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(nq)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

@qprog x_syndrome_circuit (idx) begin
    b = _xadj(idx)

    nb = length(b)
    for j in 2:nb
        CNOT(b[1], b[j])
    end
    H(b[1])
    res = M(b[1])
    H(b[1])
    for j in nb:-1:2
        CNOT(b[1], b[j])
    end

    res
end

@qprog z_syndrome_circuit (idx) begin
    b = _zadj(idx)

    nb = length(b)
    for j in 2:nb
        CNOT(b[j], b[1])
    end
    res = M(b[1])
    for j in nb:-1:2
        CNOT(b[j], b[1])
    end

    res
end

function bug(ρ, r_x, r_z, d)
    println("NO BUG....")
end

@qprog qec_decoder (nx, nz, nq, d, ctx) begin
    # println("Decoder start")

    s_x = [x_syndrome_circuit(j) for j in 1:nx]
    s_z = [z_syndrome_circuit(j) for j in 1:nz]

    r_x = decoder_algo_xz(ctx, d, s_x, "X", nq, _xadj)
    r_z = decoder_algo_xz(ctx, d, s_z, "Z", nq, _zadj)

    for j in 1:nq
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # bug(ρ, r_x, r_z, d)
    e = reduce(&, r_z[1:((d-1)÷2)])
    sX(1, e)

    # println("Decoder end")
end

# export decoder_algo_xz, qec_default_x_syn_ckt, qec_default_z_syn_ckt

end