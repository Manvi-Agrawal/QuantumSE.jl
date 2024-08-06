# How do I write decoder for no ancilla qubits?

# This fails at Symbolic Execution Stage aka "cfgs1 = QuantSymEx(cfg1)".

using QuantumSE
using Z3

ctx = Context()

_adj(idx) = [ idx, idx%2 + 1 ]

function mwpm(d::Integer, s, s_type)

    # pre-condition
    ϕ₁ = simplify(reduce(⊻, s)) == _bv_val(ctx, 0)

    # post-condition
    ϕ₂ = bool_val(ctx, true)
    
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:dq]
    for j in 1:dq
        ϕ₂ = ϕ₂ & ((s[j] ⊻ reduce(⊻, r[[_adj(j)...]])) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(2*d*d)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(2*d*d)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end


@qprog toric_x_m (idx) begin
    # b = _xadj(d, idx)
            
    # CNOT(idx, b[1])
    # CNOT(idx, b[2])

    # H(idx)
    # res = M(idx)
    # H(idx)

    # CNOT(idx, b[2])
    # CNOT(idx, b[1])

    # res
end



@qprog toric_z_m (d, idx) begin
    # b = _zadj(d, idx)
        
    # CNOT(idx, b[1])
    # CNOT(idx, b[2])

    # res = M(idx)

    # CNOT(idx, b[2])
    # CNOT(idx, b[1])

    # res
end


@qprog toric_decoder (d) begin

    # s_x = [toric_x_m(d, j) for j in 1:d*d]
    # s_z = [toric_z_m(d, j) for j in 1:d*d]
    
    # r_x = mwpm(d, s_x, "X")
    # r_z = mwpm(d, s_z, "Z")
    
    # for j in 1:d*d
    #     sZ(j, r_x[j])
    #     sX(j, r_z[j])
    # end

    # a strange bug
    # e = reduce(&, r_z[1])

    # sX(1, e)

end

function check_toric_decoder()
    # TODO: Remove hardcoding
    d = 2

    @info "Initailization Stage"
    t0 = time()
    begin
        num_qubits = d

        # 2x4 matrix
	    stabilizer = fill(false, num_qubits, 2*num_qubits)

        # X1  on D1
        stabilizer[1, 1] = true
        # Z1 on D1
        stabilizer[1, 2+1] = true


        # X2 on D2
        stabilizer[2, 2] = true
        # Z2 on D2
        stabilizer[2, 2+2] = true


	    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
	    # phases = fill(false, num_qubits)

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
            (:_adj, _adj),
            (:ctx, ctx),
            (:mwpm, mwpm)
        ])

        # num_x_errors = 0
        # x_errors = inject_errors(ρ1, "X")
        # ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

        # x_errors = inject_errors(ρ2, "X")
        # ϕ_x2 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)

        cfg1 = SymConfig(toric_decoder(d), σ, ρ1)
    end

    @info "Symbolic Execution Stage"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1)
    end

    @info "SMT Solver Stage"
    t2 = time()
    begin
        res = true
        for cfg in cfgs1
            if !check_state_equivalence(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z1=#, cfg.ϕ[1], cfg.ϕ[2]),
                `bitwuzla --smt-comp-mode true -rwl 0 -S kissat`
                #`bitwuzla --smt-comp-mode true -S kissat`
                #`bitwuzla --smt-comp-mode true -rwl 0`
               )
                res = false
                break
            end
        end
    end

    t3 = time()

    res, t3-t0, t1-t0, t2-t1, t3-t2
end


# open("toric_code.dat", "w") do io
# println(io, "nq all init qse smt")
println("nq all init qse smt")
res, all, init, qse, smt = check_toric_decoder()
println("$(nq) $(all) $(init) $(qse) $(smt)")
