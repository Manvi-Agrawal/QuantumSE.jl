using QuantumSE
using Z3

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

function toric_z_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in (_zadj(d, idx) .+ 2*d*d)
		s[j] = true
	end

	s
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

@qprog toric_decoder (d) begin

    s_x = [toric_x_m(d, j) for j in 1:d*d]
    s_z = [toric_z_m(d, j) for j in 1:d*d]
    
    r_x = mwpm(d, s_x, "X")
    r_z = mwpm(d, s_z, "Z")
    
    for j in 1:2*d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)

end

function get_stabilizer_and_phases(d::Integer)
    num_qubits = 2*d*d

    stabilizer = Matrix{Bool}(undef, num_qubits, 2*num_qubits)
    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
    lx = _bv_const(ctx, "lx")
    lz = _bv_const(ctx, "lz")

    @simd for i in 1:d*d-1
        stabilizer[i,:] = toric_x_s(d, i)
        stabilizer[i+d*d,:] = toric_z_s(d, i)
        phases[i] = _bv_val(ctx, 0)
        phases[i+d*d] = _bv_val(ctx, 0)
    end

    stabilizer[d*d,:] = toric_lx1(d)
    phases[d*d] = lx
    stabilizer[2*d*d,:] = toric_lz1(d)
    phases[2*d*d] = lz

    return (stabilizer, phases)
end

function get_config(stabilizer, phases, d)
    num_qubits = 2*d*d

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
    ρ1 = copy(ρ01)

    σ = CState([(:d, d),
        (:toric_decoder, toric_decoder),
        (:toric_z_m, toric_z_m),
        (:toric_x_m, toric_x_m),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:mwpm, mwpm)
    ])

    num_x_errors = (d-1)÷2
    x_errors = inject_errors(ρ1, "X")
    ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)


    cfg1 = SymConfig(toric_decoder(d), σ, ρ1)

    return (ρ01, ϕ_x1, cfg1)
end

function check_toric_decoder(d::Integer)

    @info "Initialization Stage"
    t0 = time()
    begin
        (stabilizer, phases) = get_stabilizer_and_phases(d)
    end

    # println("Encoded stabilizer : $(stabilizer)")

    @info "Configuration Init"
    t1 = time()
    begin
        (ρ01, ϕ_x1, cfg1) = get_config(stabilizer, phases, d)
        cfgs1 = QuantSymEx(cfg1)
    end

    @info "SMT constraint Generation"
    t2 = time()

    begin
        res = true
        for cfg in cfgs1
            if !generate_constraints(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                )
                res = false
                break
            end
        end

    end

    @info "SMT Solver Stage"
    t3 = time()

    begin
        # res = solve_constraints(`bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
    end

    t4 = time()

    res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

check_toric_decoder(3) # precompile time

open("toric_code.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 11:11#16
        res_d, all, init, config, cons_gen, cons_sol = check_toric_decoder(d)
        println("d,res,nq,all,init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(2*d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
        println(io, "$(d),$(res_d),$(2*d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
    end
end
