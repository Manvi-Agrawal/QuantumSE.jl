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




@qprog qec_decoder (ctx, ρ, d, bug) begin
    # println("Decoder start")

    nq = d*d
    lim = (nq-1)÷2

    

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
    # bug(ρ, r_x, r_z, d)

    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)



    # println("Decoder end")
end




function get_sym_config(stabilizer, decoder_config)
    d = decoder_config.d
    phases = decoder_config.phases
    X_nbr = decoder_config.X_nbr
    Z_nbr = decoder_config.Z_nbr
    ctx = decoder_config.ctx
    x_syndrome_circuit = decoder_config.x_syndrome_circuit
    z_syndrome_circuit = decoder_config.z_syndrome_circuit
    decoder_algo_xz = decoder_config.decoder_algo_xz
    bug = decoder_config.bug


    # println("LOG-INFO: Inside get_sym_config")
    num_qubits = d*d

    _xadj(j) = X_nbr[j]
    _zadj(j) = Z_nbr[j]

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)

    # println("LOG-INFO: After from_stabilizer")

    ρ1 = copy(ρ01)


    σ = CState([(:d, d),
        (:qec_decoder, qec_decoder),
        (:x_syndrome_circuit, x_syndrome_circuit),
        (:z_syndrome_circuit, z_syndrome_circuit),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:bug, bug),
        (:decoder_algo_xz, decoder_algo_xz)
    ])

    num_x_errors = (d-1)÷2
    x_errors = inject_errors(ρ1, "X")
    ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)
  
    decoder = qec_decoder(ctx, ρ1, d, bug)
    return (ρ01, ϕ_x1, SymConfig(decoder, σ, ρ1) )
end

Base.@kwdef mutable struct QecDecoderConfig
    d::Integer
    X_nbr
    Z_nbr
    stabilizer
    phases
    ctx
    x_syndrome_circuit = QEC_Defaults.x_syndrome_circuit
    z_syndrome_circuit = QEC_Defaults.z_syndrome_circuit
    decoder_algo_xz = QEC_Defaults.decoder_algo_xz
    bug = QEC_Defaults.bug
end

function check_qec_decoder(decoder_config::QecDecoderConfig)
    @info "Initialization Stage: Encode State"
    t0 = time()
    begin
        d = decoder_config.d
        X_nbr = decoder_config.X_nbr
        Z_nbr = decoder_config.Z_nbr
        stabilizer = decoder_config.stabilizer
     
        # stabilizer = get_stabilizer(d, X_nbr, Z_nbr)
        
    end

    
    # println("X nbr: $(X_nbr)")
    # println("Z nbr: $(Z_nbr)")

    # println("Encoded stabilizer : $(stabilizer)")
    # println("Phases : $(phases)")


    @info "Decoder Configuration"
    t1 = time()
    begin
        (ρ01, ϕ_x1, cfg1) = get_sym_config(stabilizer, decoder_config)
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
        res = solve_constraints()
    end

    t4 = time()

    res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

# export check_qec_decoder

end