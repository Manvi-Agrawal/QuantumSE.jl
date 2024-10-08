using QuantumSE
using Z3

N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()

function _adj4(d, idx)
    r = (idx-1)÷d
    c = (idx-1)%d

    return [r*d+c, r*d+c+1, (r+1)*d+c, (r+1)*d+(c+1)] .+ 1
end

function _adj2_hor(d, idx)
    r = (idx-1)÷d
    c = (idx-1)%d

    return [r*d + c, r*d + c+1 ] .+ 1
end

function _adj2_ver(d, idx)
    r = (idx-1)÷d
    c = (idx-1)%d

    return [r*d + c, (r+1)*d + c] .+ 1
end

function css_check(d, s, s_type, nq, adj)

    ## pre-condition
    ϕ₁ = bool_val(ctx, true)

    ## post-condition
    ϕ₂ = bool_val(ctx, true)
    r = [_bv_const(ctx, "r_$(s_type)_$(j)") for j in 1:nq]
    for j in 1:length(s)
        ϕ₂ = ϕ₂ & (s[j] ⊻ reduce(⊻,  r[adj(j)]) == _bv_val(ctx, 0))
    end

    ϕ₃ = (sum( (x -> concat(bv_val(ctx, 0, _len2(nq)), x)).(r) ) <= bv_val(ctx, (d-1)÷2, _len2(nq)+1))

    (r, ϕ₁, ϕ₂ & ϕ₃)
end


@qprog surface_code_x_m (idx) begin
    b = _xadj(idx)

    nb = length(b)

    for j in 2:nb
        CNOT(b[1], b[j])
    end

    H(b[1])
    res = M(b[1])
    H(b[1])
    
    for j in nb:-1:2
        CNOT(b[1], b[j])
    end

    res
end

@qprog surface_code_z_m (idx) begin
    b = _zadj(idx)

    nb = length(b)
    for j in 2:nb
        CNOT(b[j], b[1])
    end
    res = M(b[1])
    for j in nb:-1:2
        CNOT(b[j], b[1])
    end

    res
end


@qprog surface_code_decoder (d) begin
    # print("Decoder start")

    nq = d*d
    lim = (nq-1)÷2

    s_x = [surface_code_x_m(j) for j in 1:lim]
    s_z = [surface_code_z_m(j) for j in 1:lim]

    r_x = css_check(d, s_x, "X", nq, _xadj)
    r_z = css_check(d, s_z, "Z", nq, _zadj)

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

function get_nbr(d::Integer)
    num_qubits = d*d
    gen_s = (num_qubits-1)÷2

    X_nbr = [ [] for _ in 1:gen_s]
    Z_nbr = [ [] for _ in 1:gen_s]

    z_4q = [ 1 + r*d+c+r%2 for r in 0:(d-2) for c in 0:2:(d-2) ]

    x_4q = [ 1 + r*d+c+(r+1)%2 for r in 0:(d-2) for c in 0:2:(d-2) ]


    x_2q_up = [ 1+c for c in 0:2:(d-2) ]
    x_2q_down = [ 1+ (d-1)*d + c for c in 1:2:(d-2) ]
    x_2q = vcat(x_2q_up, x_2q_down)



    z_2q_right = [ 1 + r*d+ (d-1) for r in 0:2:(d-2) ]
    z_2q_left = [ 1 + r*d for r in 1:2:(d-2) ]
    z_2q = vcat(z_2q_left, z_2q_right)


    x_q = vcat(x_4q, x_2q)
    z_q = vcat(z_4q, z_2q)


    for j in eachindex(x_4q)
        for x_adj in  _adj4( d, x_4q[j])
            push!(X_nbr[j], x_adj)
        end
    end

    for j in eachindex(x_2q)
        for x_adj in _adj2_hor( d, x_2q[j])
            push!(X_nbr[length(x_4q) + j], x_adj)
        end
    end

    # r += length(x_2q)

    for j in eachindex(z_4q)
        for z_adj in _adj4( d, z_4q[j])
            push!(Z_nbr[j], z_adj)
        end
    end

    for j in eachindex(z_2q)
        for z_adj in _adj2_ver( d, z_2q[j])
            push!(Z_nbr[length(z_4q) + j], z_adj)
        end
    end

    

    return (X_nbr, Z_nbr)

end

function get_config(stabilizer, phases, X_nbr, Z_nbr, d::Integer)
    num_qubits = d*d

    _xadj(j) = X_nbr[j]
    _zadj(j) = Z_nbr[j]

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
    ρ1 = copy(ρ01)

    σ = CState([(:d, d),
        (:surface_code_decoder, surface_code_decoder),
        (:surface_code_z_m, surface_code_z_m),
        (:surface_code_x_m, surface_code_x_m),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:css_check, css_check)
    ])

    num_x_errors = (d-1)÷2
    x_errors = inject_errors(ρ1, "X")
    ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)
  
    decoder = surface_code_decoder(d)
    return (ρ01, ϕ_x1, SymConfig(decoder, σ, ρ1) )
end

function get_phases(d::Integer)
    num_qubits = d*d

    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
    lx = _bv_const(ctx, "lx")

    for i in 1:num_qubits-1
        phases[i] = _bv_val(ctx, 0)
    end

    phases[num_qubits] = lx

    return phases
end

function check_surface_code_decoder(d::Integer)
    @info "Initialization Stage: Encode State"
    t0 = time()
    begin
        ( X_nbr, Z_nbr) = get_nbr(d)
        stabilizer = get_stabilizer(d, X_nbr, Z_nbr)
        phases = get_phases(d)
        (ρ01, ϕ_x1, cfg1) = get_config(stabilizer, phases, X_nbr, Z_nbr, d)
    end

    
    # println("X idxs: $(X_nbr)")
    # println("Z idxs: $(Z_nbr)")

    # println("Encoded stabilizer : $(stabilizer)")

    @info "Decoder Configuration"
    t1 = time()
    begin
        cfgs1 = QuantSymEx(cfg1)
    end

    @info "Constraint Generation"
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

    @info "Constraint Solver"
    t3 = time()

    begin
        res = solve_constraints(`bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
    end

    t4 = time()

    res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

check_surface_code_decoder(3) # precompile time
@info "precompile done..."

open("surface_code.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 3:2:11
        res_d, all, init, config, cons_gen, cons_sol = check_surface_code_decoder(d)
        println("d,res,nq,all,init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
        println(io, "$(d),$(res_d),$(d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
    end
end
