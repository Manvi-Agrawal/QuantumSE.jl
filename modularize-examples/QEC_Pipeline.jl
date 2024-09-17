module QEC_Pipeline



using QuantumSE
using Z3

include("QEC_Defaults.jl")
using .QEC_Defaults

include("QEC_Helper.jl")
using .QEC_Helper


N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()

function get_sym_config(stabilizer, decoder_config)
    d = decoder_config.d
    phases = decoder_config.phases
    _xadj = decoder_config._xadj
    _zadj = decoder_config._zadj
    ctx = decoder_config.ctx
    x_syndrome_circuit = decoder_config.x_syndrome_circuit
    z_syndrome_circuit = decoder_config.z_syndrome_circuit
    decoder_algo_xz = decoder_config.decoder_algo_xz
    bug = decoder_config.bug


    # println("LOG-INFO: Inside get_sym_config")
    num_qubits = decoder_config.num_qubits

    

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)

    # println("LOG-INFO: After from_stabilizer")

    ρ1 = copy(ρ01)

    qec_decoder = decoder_config.decoder
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
  
    nq = num_qubits
    xlim = (nq-1)÷2
    zlim = (nq-1)÷2

    decoder = qec_decoder(decoder_config.decoder_params...)
    return (ρ01, ϕ_x1, SymConfig(decoder, σ, ρ1) )
end

Base.@kwdef mutable struct QecPipelineConfig
    d::Integer
    num_qubits::Integer
    stabilizer
    phases
    ctx
    _xadj = missing
    _zadj = missing
    decoder = QEC_Defaults.qec_decoder
    decoder_params = missing
    x_syndrome_circuit = QEC_Defaults.x_syndrome_circuit
    z_syndrome_circuit = QEC_Defaults.z_syndrome_circuit
    decoder_algo_xz = QEC_Defaults.decoder_algo_xz
    bug = QEC_Defaults.bug
end

function check_qec_decoder(decoder_config::QecPipelineConfig)
    @info "Initialization Stage"
    t0 = time()
    begin
        d = decoder_config.d
        nq = decoder_config.num_qubits
        stabilizer = decoder_config.stabilizer
        (ρ01, ϕ_x1, cfg1) = get_sym_config(stabilizer, decoder_config)
    end
    

    # Res should be true initially, set false to disable later code
    res = true

    begin
        t1 = time()
        @info "Symbolic Execution"
        cfgs1 = res ? QuantSymEx(cfg1) : [cfg1]
    end

    begin
        t2 = time()
        @info "Constraint Generation"

        for cfg in cfgs1
            if res && !generate_constraints(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                )
                res = false
                break
            end
        end
    end

    begin
        t3 = time()
        @info "Constraint Solver"
        res = res && solve_constraints()
    end

    t4 = time()

    res, d, nq, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

function qec_runner(file, get_decoder_config, range)
    qec_pre = get_decoder_config(range[1])
    QEC_Pipeline.check_qec_decoder(qec_pre)
    @info "PRECOMPLIED qec code..."

    open(file, "w") do io
        println(io, "d,res,nq,all,init_config,init,config,cons_gen,cons_sol")

    
        for r in range
            tm2 = time()
    
            decoder_config = get_decoder_config(r...)
    
            tm1 = time()
    
            res_d, d, nq, all, init, config, cons_gen, cons_sol = QEC_Pipeline.check_qec_decoder(decoder_config)
    
            init_config = (tm1-tm2)
            all += init_config
    
            println("d,res,nq,all,init_config,init,config,cons_gen,cons_sol")
            println("$(d),$(res_d),$(nq),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
            println(io, "$(d),$(res_d),$(nq),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
        end
    end    
end

end