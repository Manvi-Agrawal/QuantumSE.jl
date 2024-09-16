module QEC_Pipeline



using QuantumSE
using Z3

include("QEC_Defaults.jl")
using .QEC_Defaults


N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()



function get_stabilizer(d::Integer, X_nbr, Z_nbr)
    num_qubits = d*d

    stabilizer = fill(false, num_qubits, 2*num_qubits)

    for row in 1:length(X_nbr)
        for col in 1:length(X_nbr[row])
            s_col = X_nbr[row][col]
            # println("Set $(row), $(s_col) to true...")
            stabilizer[row, s_col] = true
        end
    end
    
    for row in 1 : length(Z_nbr)
        for col in 1:length(Z_nbr[row])
            s_row = row + length(X_nbr)
            s_col = num_qubits+Z_nbr[row][col]
            stabilizer[s_row, s_col] = true
        end
    end 

    # r: 1, 2, ..., d ;; c = 1
    sc_lx = [ 1+ r*d + 1 for r in 0:(d-1) ]

    for idx in sc_lx
        stabilizer[num_qubits, idx] = true
    end

    return stabilizer

end

@qprog qec_decoder (ctx, d) begin
    # print("Decoder start")

    nq = d*d
    lim = (nq-1)÷2

    x_syndrome_circuit = qec_default_x_syn_ckt
    z_syndrome_circuit = qec_default_z_syn_ckt

    s_x = [x_syndrome_circuit(j) for j in 1:lim]
    s_z = [z_syndrome_circuit(j) for j in 1:lim]

    r_x = decoder_algo_xz(ctx, d, s_x, "X", nq, _xadj)
    r_z = decoder_algo_xz(ctx, d, s_z, "Z", nq, _zadj)

    # print("rx=$(r_x)")
    # print("rz=$(r_z)")


    for j in 1:d*d
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # e = reduce(&, r_z[1:(d-1)÷2])

    # sX(1, e)

    # print("Decoder end")
end




function get_config(stabilizer, phases, X_nbr, Z_nbr, d::Integer, ctx)
    println("LOG-INFO: Inside get_config")
    num_qubits = d*d

    _xadj(j) = X_nbr[j]
    _zadj(j) = Z_nbr[j]

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)

    println("LOG-INFO: After from_stabilizer")

    ρ1 = copy(ρ01)

    σ = CState([(:d, d),
        (:qec_decoder, qec_decoder),
        (:qec_default_x_syn_ckt, qec_default_x_syn_ckt),
        (:qec_default_z_syn_ckt, qec_default_z_syn_ckt),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:decoder_algo_xz, decoder_algo_xz)
    ])

    num_x_errors = (d-1)÷2
    x_errors = inject_errors(ρ1, "X")
    ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)
  
    decoder = qec_decoder(ctx, d)
    return (ρ01, ϕ_x1, SymConfig(decoder, σ, ρ1) )
end



function check_qec_decoder(d::Integer, X_nbr, Z_nbr, phases, ctx)
    @info "Initialization Stage: Encode State"
    t0 = time()
    begin
        stabilizer = get_stabilizer(d, X_nbr, Z_nbr)
    end

    
    # println("X nbr: $(X_nbr)")
    # println("Z nbr: $(Z_nbr)")

    println("Encoded stabilizer : $(stabilizer)")
    println("Phases : $(phases)")


    @info "Decoder Configuration"
    t1 = time()
    begin
        (ρ01, ϕ_x1, cfg1) = get_config(stabilizer, phases, X_nbr, Z_nbr, d, ctx)
        cfgs1 = QuantSymEx(cfg1)
    end

    @info "Constraint Generation"
    t2 = time()

    begin
        res = true
        for cfg in cfgs1
            if res && !generate_constraints(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                )
                res = false
                break
            end
        end

    end

    @info "Constraint Solver"
    t3 = time()

    begin
        res = solve_constraints(`bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
    end

    t4 = time()

    res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

export check_qec_decoder

end