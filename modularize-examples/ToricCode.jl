using QuantumSE
using Z3

# include("QEC_Helper.jl")
# using .QEC_Helper

include("QEC_Pipeline.jl")
using .QEC_Pipeline

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

    println("Type(phases): $(typeof(phases))")
    println("Length(phases): $(length(phases))")

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

@qprog toric_decoder (d) begin

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

# function toric_decoder(d, ctx, _xadj, _zadj)
#     nq = 2*d*d
#     xlim = (nq-1)÷2
#     zlim = (nq-1)÷2

#     decoder = QEC_Helper.qec_decoder(ctx, d, nq, xlim, zlim, _xadj, _zadj, bug)
# end


open("toric_code.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 3:7
        tm2 = time()
        # (X_nbr, Z_nbr) = get_nbr(d)

        stabilizer = get_stabilizer(d)

        toric_decoder_params = (d)
        toric_decoder_config = QEC_Pipeline.QecDecoderConfig(
            d=d,
            num_qubits=2*d*d,
            stabilizer=stabilizer,
            phases= get_phases(d),
            ctx=ctx,
            _xadj=_xadj,
            _zadj=_zadj,
            decoder=toric_decoder,
            decoder_params=toric_decoder_params)

        toric_decoder_config.x_syndrome_circuit = toric_x_m
        toric_decoder_config.z_syndrome_circuit = toric_z_m
        toric_decoder_config.decoder_algo_xz = mwpm

        


        tm1 = time()
        

        # println("X_nbr: $(X_nbr)")
        # println("Z_nbr: $(Z_nbr)")

        # println("Phases: $(phases)")

        res_d, all, init, config, cons_gen, cons_sol = QEC_Pipeline.check_qec_decoder(toric_decoder_config)

        init_config = (tm2-tm1)
        all += init_config

        println("d,res,nq,all,init_config, init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(2*d*d),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
        println(io, "$(d),$(res_d),$(2*d*d),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
    end
end
