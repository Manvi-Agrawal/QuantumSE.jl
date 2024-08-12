# This example fails on from_stabilizer call
# ERROR: LoadError: DimensionMismatch: number of columns of each array must match (got (8, 4))

using QuantumSE
using Z3

ctx = Context()

function _xadj(d, q_idx)
    if q_idx == 1
        return [2]
    elseif q_idx == 2
        return [1, 4]
    elseif q_idx == 4
        return [2]
    else
        return [3]
    end
end

function _zadj(d, q_idx)
    if q_idx == 1
        return [3]
    elseif q_idx == 3
        return [1, 4]
    elseif q_idx == 4
        return [3]
    else
        return []
    end
end



@qprog toric_decoder (d) begin

    # aq_x = [2]
    # s_x = [toric_x_m(d, j) for j in aq_x]

    # aq_z = [3]
    # s_z = [toric_z_m(d, j) for j in aq_z]
    
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

        # # X1  on D1
        # stabilizer[1, 0+1] = true
        # # Z1 on D1
        # stabilizer[1, 2+1] = true

        # # X1  on A1(aka q2)
        # stabilizer[2, 0+1] = true
        # # X2 on A1
        # stabilizer[2, 0+2] = true

        # # Z1  on A2(aka q3)
        # stabilizer[3, 2+1] = true
        # # Z2 on A2
        # stabilizer[3, 2+2] = true

        # # X2 on D2
        # stabilizer[4, 0+2] = true
        # # Z2 on D2
        # stabilizer[4, 2+2] = true

	    # phases = Vector{Z3.ExprAllocated}(undef, num_qubits)

        # phases[i] = _bv_val(ctx, 0)
	    phases = fill( _bv_val(ctx, 0), num_qubits)

	    lx = _bv_const(ctx, "lx")
	    lz = _bv_const(ctx, "lz")

        print("Stabilizer init: $(stabilizer)")

        # [
        # X1 X2 Z1 Z2
        # 1  0  1  0; D1
        # 1  1  0  0; A1
        # 0  0  1  1; A2
        # 0  1  0  1; D2
        # ]

         # [
        # X1 X2 X3 X4     Z1 Z2 Z3 Z4
        # 1  1  1  1      0  0  0  0  
        # 0  0  0  0      1  1  1  1  
        # 1  0  1  0      0  0  0  0  
        # 0  0  0  0      0  1  0  1  
        # ]

        
        # [
        # X1 X2 X3 X4  Z1 Z2 Z3 Z4
        # 1  0  0 0     1  0 0 0   D1
        # 1  1  0 0     0  0 0 0   A1
        # 0  0  0 0     1  1 0 0   A2
        # 0  1  0 0     0  1 0 0   D2
        # ]

        ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
        ρ1 = copy(ρ01)



        σ = CState([(:d, d),
            (:toric_decoder, toric_decoder),
            # (:toric_z_m, toric_z_m),
            # (:toric_x_m, toric_x_m),
            # (:_xadj, _xadj),
            # (:_zadj, _zadj),
            (:ctx, ctx),
            # (:mwpm, mwpm)
        ])

        # σ = CState([
        # ])

        num_x_errors = 0
        x_errors = inject_errors(ρ1, "X")
        ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)



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
d=2
res, all, init, qse, smt = check_toric_decoder(d)
println("d all init qse smt")
println("$(d) $(all) $(init) $(qse) $(smt)")
