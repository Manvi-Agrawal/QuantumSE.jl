# This example fails on from_stabilizer call
# ERROR: LoadError: DimensionMismatch: number of columns of each array must match (got (8, 4))

using QuantumSE
using Z3

ctx = Context()

function _xadj(d, q_idx)
    if q_idx == 2
        return [1, 4]
    else
        return [3]
    end
end

function _zadj(d, q_idx)
    if q_idx == 3
        return [1, 4]
    else
        return []
    end
end

function mwpm(d::Integer, s, s_type)

    dq = [1, 4]

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    adj = s_type == "X" ? _xadj : _zadj
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in dq]
    for j in dq
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[adj(d, j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end



@qprog toric_x_m (d, idx) begin
    b = _xadj(d, idx)
    
    H(idx)
    CNOT(idx, b[1])
    CNOT(idx, b[2])
    H(idx)

    res = M(idx)

    H(idx)
    CNOT(idx, b[2])
    CNOT(idx, b[1])
    H(idx)

    
    res
end



@qprog toric_z_m (d, idx) begin
    b = _zadj(d, idx)
    
    CNOT(b[1], idx)
    CNOT(b[2], idx)

    res = M(idx)

    CNOT(b[2], idx)
    CNOT(b[1], idx)

    res
end


@qprog toric_decoder (d) begin

    aq_x = [2]
    s_x = [toric_x_m(d, j) for j in aq_x]

    aq_z = [3]
    s_z = [toric_z_m(d, j) for j in aq_z]
    
    r_x = mwpm(d, s_x, "X")
    r_z = mwpm(d, s_z, "Z")
    
    for j in 1:d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # e = reduce(&, r_z[1])

    # sX(1, e)

end

function check_toric_decoder(d::Integer)
    # TODO: Remove hardcoding
    d = 2

    @info "Initailization Stage"
    t0 = time()
    begin
        num_qubits = d*d

        # 4x4 matrix
	    # stabilizer = fill(false, num_qubits, 2*d)

        stabilizer = [
        true   true   true   true      false  false  false  false;  
        false  false  false  false     true   true   true   true;  
        true   false  true  false      false  false  false  false;  
        false  false  false  false      false  true  false  true]

	    phases = fill( _bv_val(ctx, 0), num_qubits)

	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

        print("Stabilizer init: $(stabilizer)")

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

        # σ = CState([
        # ])

        num_x_errors = 0
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)



        # x_errors = inject_errors(ρ2, "X")
        # ϕ_x2 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

        decoder = toric_decoder(d)

        println("Decoder created")

        cfg1 = SymConfig(decoder, σ, ρ1)

        println("Sym config created")

    end

    @info "Symbolic Execution Stage"
    t1 = time()
    # begin
    #     cfgs1 = QuantSymEx(cfg1)
    # end

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

    res, t3-t0, t1-t0, t2-t1, t3-t2
end


# open("toric_code.dat", "w") do io
# println(io, "nq all init qse smt")
d=2
res, all, init, qse, smt = check_toric_decoder(d)
println("d all init qse smt")
println("$(d) $(all) $(init) $(qse) $(smt)")
