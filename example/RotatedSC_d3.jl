using QuantumSE
using Z3

ctx = Context()

function _xadj(idx, d)
    if idx == 1
        return [1, 2]
    elseif idx == 2
            return [2, 3, 5, 6]
    elseif idx == 4
        return [4, 5, 7, 8]

    elseif idx == 8
        return [8, 9]
    else
        return [] 
    end
    
end

function _zadj(idx, d)
    if idx == 1
        return [1, 2, 4, 5]
    elseif idx == 3
            return [3, 6]
    elseif idx == 4
        return [4, 7]
    elseif idx == 5
        return [5, 6, 8, 9]
    else
        return [] 
    end
    
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

@qprog surface_code_z_m (d, idx) begin
    b = _zadj(d, idx)

    if length(b) == 4
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[2], b[4])
        res = M(b[4])
        CNOT(b[2], b[4])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])
    elseif length(b) == 2
        CNOT(b[1], b[2])
        res = M(b[1])
        CNOT(b[1], b[2])
    else
        res = -100
    end
    
    res
end


@qprog surface_code_x_m (d, idx) begin
    b = _xadj(d, idx)

    if length(b) == 4
        CNOT(b[1], b[2])
        CNOT(b[3], b[4])
        CNOT(b[3], b[1])
        H(b[3])
        res = M(b[3])
        H(b[3])
        CNOT(b[3], b[1])
        CNOT(b[3], b[4])
        CNOT(b[1], b[2])

    elseif length(b) == 2
        CNOT(b[1], b[2])
        H(b[1])
        res = M(b[1])
        H(b[1])
        CNOT(b[1], b[2])
    else
        res = -100
    end
    
    res
end

@qprog surface_code_decoder (d) begin

    xq = [1, 2, 4, 8]
    zq = [1, 3, 4, 5]

    s_x = [surface_code_x_m(d, j) for j in xq]
    s_z = [surface_code_z_m(d, j) for j in zq]
    
    r_x = mwpm(d, s_x, "X")
    r_z = mwpm(d, s_z, "Z")
    
    for j in 1:2*d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # e = reduce(&, r_z[1:(d-1)÷2])

    # sX(1, e)

end

function toric_x_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in _xadj(d, idx)
		s[j] = true
	end

	s
end

function check_surface_code_decoder(d::Integer)

    d = 3

    @info "Initailization Stage"
    t0 = time()
    begin
        num_qubits = d*d*2

	    stabilizer = fill(false, num_qubits, 2*num_qubits)

        xq = [1, 2, 4, 8]

        for r in 1:4
            for x_adj in _xadj(xq[r], 3)
                stabilizer[r, x_adj] = true
            end
        end

        zq = [1, 3, 4, 5]

        for r in 1:4
            for z_adj in _zadj(zq[r], 3)
                stabilizer[r, z_adj] = true
            end
        end



	    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

	    @simd for i in 1:d*d-1
	    	# stabilizer[i,:] = toric_x_s(d, i)
	    	# stabilizer[i+d*d,:] = toric_z_s(d, i)
	    	phases[i] = _bv_val(ctx, 0)
	    	phases[i+d*d] = _bv_val(ctx, 0)
	    end

	    # stabilizer[d*d,:] = toric_lx1(d)
	    phases[d*d] = lx
	    # stabilizer[2*d*d,:] = toric_lz1(d)
	    phases[2*d*d] = lz

        # print("Encoded Stabilizer 1: $(stabilizer)\n")

        ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ1 = copy(ρ01)

        # stabilizer[d*d,:] = toric_lx2(d)
	    # phases[d*d] = lx
	    # stabilizer[2*d*d,:] = toric_lz2(d)
	    # phases[2*d*d] = lz

        # print("Encoded Stabilizer 2: $(stabilizer)\n")

        # ρ02 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        # ρ2 = copy(ρ02)

        σ = CState([(:d, d),
            (:surface_code_decoder, surface_code_decoder),
            (:surface_code_z_m, surface_code_z_m),
            (:surface_code_x_m, surface_code_x_m),
            (:_xadj, _xadj),
            (:_zadj, _zadj),
            (:ctx, ctx),
            (:mwpm, mwpm)
        ])

        num_x_errors = (d-1)÷2
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

	  
        decoder = surface_code_decoder(d)
        cfg1 = SymConfig(decoder, σ, ρ1)
        # cfg2 = SymConfig(surface_code_decoder(d), σ, ρ2)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1)
        # cfgs2 = QuantSymEx(cfg2)
    end

    @info "SMT Solver Stage"
    t2 = time()

    begin
        res = true
        # for cfg in cfgs1
        if !check_state_equivalence(
            ρ01, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg1.ϕ[1], cfg1.ϕ[2]),
            `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
            #`bitwuzla --smt-comp-mode true -S kissat`
            #`bitwuzla --smt-comp-mode true -rwl 0`
            )
            res = false
            # break
        end
        # end

    end

    t3 = time()

    res, t2-t0, t1-t0, t2-t1, t3-t2
end

# check_surface_code_decoder(3) # precompile time

d = 3
res, all, init, qse, smt = check_surface_code_decoder(d)
println("d: res nq all init qse smt")
println("$(d): $(res) $(d*d*2) $(all) $(init) $(qse) $(smt)")

