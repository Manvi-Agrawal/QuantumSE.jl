using QuantumSE
using Z3

include("QEC_Pipeline.jl")
using .QEC_Pipeline

include("QEC_Decoder_Defaults.jl")
using .QEC_Decoder_Defaults

include("QecDecoder.jl")
include("QecPipelineConfig.jl")

ctx = Context()

function _xadj(d, idx)
    ai = (idx-1)÷d
    bi = (idx-1)%d

    [d*d+ai*d+bi, ai*d+(bi+d-1)%d, ai*d+bi, d*d+((ai+d-1)%d)*d+bi] .+ 1
end

function _zadj(d, idx)
    ai = (idx-1)÷d
    bi = (idx-1)%d

    [((ai+1)%d)*d+bi, d*d+ai*d+bi, d*d+ai*d+(bi+1)%d, ai*d+bi] .+ 1
end

function mwpm(d::Integer, s, s_type)

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:2*d*d]
    for j in 1:d*d
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[adj(d, j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end

function toric_x_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in _xadj(d, idx)
		s[j] = true
	end

	s
end

function toric_z_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in (_zadj(d, idx) .+ 2*d*d)
		s[j] = true
	end

	s
end

function toric_lx1(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[i*d] = true
	end

	s
end

function toric_lx2(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[d*d+i] = true
	end

	s
end

function toric_lz1(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[3*d*d+i*d] = true
	end

	s
end

function toric_lz2(d::Integer)
	s = zeros(Bool, 4*d*d)

	@inbounds @simd for i in 1:d
		s[2*d*d+i] = true
	end

	s
end

@qprog toric_x_m (d, idx) begin
    b = _xadj(d, idx)
    
    CNOT(b[1], b[2])
    CNOT(b[3], b[4])
    CNOT(b[3], b[1])
    H(b[3])
    res = M(b[3])
    H(b[3])
    CNOT(b[3], b[1])
    CNOT(b[3], b[4])
    CNOT(b[1], b[2])
    
    res
end



@qprog toric_z_m (d, idx) begin
    b = _zadj(d, idx)
    
    CNOT(b[1], b[2])
    CNOT(b[3], b[4])
    CNOT(b[2], b[4])
    res = M(b[4])
    CNOT(b[2], b[4])
    CNOT(b[3], b[4])
    CNOT(b[1], b[2])
    
    res
end



function get_phases(d)
    num_qubits = 2*d*d

    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
    lx = _bv_const(ctx, "lx")
    lz = _bv_const(ctx, "lz")
    
    for i in 1:d*d-1
        phases[i] = _bv_val(ctx, 0)
        phases[i+d*d] = _bv_val(ctx, 0)
    end

    phases[d*d] = lx
    phases[2*d*d] = lz

    return phases
end

function get_stabilizer(d)
    num_qubits = 2*d*d
    stabilizer = Matrix{Bool}(undef, num_qubits, 2*num_qubits)

    for i in 1:d*d-1
        stabilizer[i,:] = toric_x_s(d, i)
        stabilizer[i+d*d,:] = toric_z_s(d, i)
    end

    stabilizer[d*d,:] = toric_lx1(d)
    stabilizer[2*d*d,:] = toric_lz1(d)

    println("Type(stabilizer): $(typeof(stabilizer))")

    return stabilizer
end

@qprog toric_decoder_ckt (d) begin

    s_x = [x_syndrome_circuit(d, j) for j in 1:d*d]
    s_z = [z_syndrome_circuit(d, j) for j in 1:d*d]
    
    r_x = decoder_algo_xz(d, s_x, "X")
    r_z = decoder_algo_xz(d, s_z, "Z")
    
    for j in 1:2*d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)

end

function get_toric_decoder(d, ctx)
    return QecDecoder(
        d=d,
        num_qubits=2*d*d,
        ctx=ctx,
        _xadj=_xadj,
        _zadj=_zadj,
        nx=d*d,
        nz=d*d,
        decoder = toric_decoder_ckt,
        decoder_params = (d),
        x_syndrome_circuit = toric_x_m,
        z_syndrome_circuit = toric_z_m,
        decoder_algo_xz = mwpm
    )

end

function get_toric_decoder_config(d)
    stabilizer = get_stabilizer(d)


    toric_decoder_config = QecPipelineConfig(
        stabilizer=stabilizer,
        phases= get_phases(d),
        decoder=get_toric_decoder(d, ctx)
        )

    # toric_decoder_config.x_syndrome_circuit = toric_x_m
    # toric_decoder_config.z_syndrome_circuit = toric_z_m
    # toric_decoder_config.decoder_algo_xz = mwpm
    # toric_decoder_config.decoder = toric_decoder
    # toric_decoder_config.decoder_params = (d)

    return toric_decoder_config
end

QEC_Pipeline.qec_runner("toric_code.csv", get_toric_decoder_config, 3:10)
