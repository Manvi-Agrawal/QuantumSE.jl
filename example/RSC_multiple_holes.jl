using QuantumSE
using Z3

using SimpleGF2
using LinearAlgebra
using InvertedIndices

function delete_matrix_cols(a, cols)
    a[:, Not(cols)]
end

function delete_row_cols(a, cols)
    a[:, Not(cols)]
end

N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()



function _adj_hole(R, C, idx)
    r = (idx-1)÷C
    c = (idx-1)%C

    valid = (idx, R, C) -> (idx>0 && idx<R*C)

    N = (r-1)*C+c + 1
    S = (r+1)*C+c + 1
    W = r*C+(c-1) + 1
    E = r*C+(c+1) + 1

    NW = (r-1)*C+(c-1) + 1
    NE = (r-1)*C+(c+1) + 1
    SW = (r+1)*C+(c-1) + 1
    SE = (r+1)*C+(c+1) + 1

    adj = [N, S, W, E, NW, NE, SW, SE]
    # print("adj: $(adj)")


    return [ x for x in adj if valid(x, R, C)] 
end

function _adj4(C, idx)
    r = (idx-1)÷C
    c = (idx-1)%C

    return [r*C+c, r*C+c+1, (r+1)*C+c, (r+1)*C+(c+1)] .+ 1
end

function _adj2_hor(C, idx)
    r = (idx-1)÷C
    c = (idx-1)%C

    return [r*C + c, r*C + c+1 ] .+ 1
end

function _adj2_ver(C, idx)
    r = (idx-1)÷C
    c = (idx-1)%C

    return [r*C + c, (r+1)*C + c] .+ 1
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


@qprog surface_code_decoder (R, C) begin
    # print("Decoder start")

    nq = R*C
    x_gen_s = ((R-1)*(C-1))÷2 + (C-1)
    z_gen_s = ((R-1)*(C-1))÷2 + (R-1)

    s_x = [surface_code_x_m(j) for j in 1:x_gen_s]
    s_z = [surface_code_z_m(j) for j in 1:z_gen_s]

    d = min(R, C)

    r_x = css_check(d, s_x, "X", nq, _xadj)
    r_z = css_check(d, s_z, "Z", nq, _zadj)

    # print("rx=$(r_x)")
    # print("rz=$(r_z)")


    for j in 1:nq
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    e = reduce(&, r_z[1:(d-1)÷2])

    sX(1, e)

    # print("Decoder end")
end



function get_nbr(R::Integer, C::Integer)
    num_qubits = R*C

    x_gen_s = ((R-1)*(C-1))÷2 + (C-1)
    z_gen_s = ((R-1)*(C-1))÷2 + (R-1)


    stabilizer = fill(false, num_qubits, 2*num_qubits)

    X_nbr = [ [] for _ in 1:x_gen_s]
    Z_nbr = [ [] for _ in 1:z_gen_s]

    z_4q = [ 1 + r*C+c+r%2 for r in 0:(R-2) for c in 0:2:(C-2) ]
    # println("Z 4q: $(z_4q)")

    x_4q = [ 1 + r*C+c+(r+1)%2 for r in 0:(R-2) for c in 0:2:(C-2) ]
    # println("X 4q: $(x_4q)")


    x_2q_up = [ 1+c for c in 0:2:(C-2) ]
    x_2q_down = [ 1+ (R-1)*C + c for c in 1:2:(C-2) ]
    x_2q = vcat(x_2q_up, x_2q_down)

#    println("X 2q: $(x_2q)")


    z_2q_right = [ 1 + r*C+ (C-1) for r in 0:2:(R-2) ]
    z_2q_left = [ 1 + r*C for r in 1:2:(R-2) ]
    z_2q = vcat(z_2q_left, z_2q_right)

#    println("Z 2q: $(z_2q)")

    x_q = vcat(x_4q, x_2q)
    z_q = vcat(z_4q, z_2q)



    r = 0

    for j in eachindex(x_4q)
        for x_adj in  _adj4( C, x_4q[j])
            stabilizer[r+j, x_adj] = true
            push!(X_nbr[j], x_adj)
        end
    end

    r += length(x_4q)


    for j in eachindex(x_2q)
        for x_adj in _adj2_hor( C, x_2q[j])
            stabilizer[r+j, x_adj] = true
            # println("Push X_nbr[$(length(x_4q) + j)] <- $(x_adj)")

            push!(X_nbr[length(x_4q) + j], x_adj)
        end
    end

    r += length(x_2q)

    for j in eachindex(z_4q)
        for z_adj in _adj4( C, z_4q[j])
            stabilizer[r+j, num_qubits+z_adj] = true
            push!(Z_nbr[j], z_adj)
        end
    end

    r += length(z_4q)

    for j in eachindex(z_2q)
        for z_adj in _adj2_ver( R, z_2q[j])
            stabilizer[r+j, num_qubits+z_adj] = true
            push!(Z_nbr[length(z_4q) + j], z_adj)
        end
    end

    return ( X_nbr, Z_nbr)

end

function get_config(stabilizer, phases, X_nbr, Z_nbr, R::Integer, C::Integer)
    num_qubits = R*C

    _xadj(j) = X_nbr[j]
    _zadj(j) = Z_nbr[j]

    ρ01 = from_stabilizer(num_qubits, stabilizer, phases, ctx)
    ρ1 = copy(ρ01)

    σ = CState([(:R, R),
        (:C, C),
        (:surface_code_decoder, surface_code_decoder),
        (:surface_code_z_m, surface_code_z_m),
        (:surface_code_x_m, surface_code_x_m),
        (:_xadj, _xadj),
        (:_zadj, _zadj),
        (:ctx, ctx),
        (:css_check, css_check)
    ])

    d = min(R, C)

    num_x_errors = (d-1)÷2
    x_errors = inject_errors(ρ1, "X")
    ϕ_x1 = _sum(ctx, x_errors, num_qubits) == bv_val(ctx, num_x_errors, _len2(num_qubits)+1)
  
    decoder = surface_code_decoder(R, C)
    return (ρ01, ϕ_x1, SymConfig(decoder, σ, ρ1) )
end

function get_phases(R::Integer, C::Integer)
    num_qubits = R*C

    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)
    lx = _bv_const(ctx, "lx")

    for i in 1:num_qubits-1
        phases[i] = _bv_val(ctx, 0)
    end

    phases[num_qubits] = lx

    return phases
end

function get_logical_op_holes(R, C, holes, l_op, lx_type)
    nq = R*C

    s_col_prefix = (lx_type=="Z") ? nq : 0

    stabilizer_holes = fill(false, length(holes), 2*nq)

    for row in 1:length(holes)
        nbr = _adj_hole(R, C, holes[row])
        println("Hole nbr: $(nbr)")

        for col in 1:length(nbr)
            s_col = s_col_prefix + nbr[col]
            stabilizer_holes[row, s_col] = true
            l_op[s_col] = false
        end
    end

    println("stabilizer_holes: $(stabilizer_holes)")
    println("l_op: $(l_op)")

    temp1 = vcat(stabilizer_holes, l_op')

    del_col = vcat(holes, [nq+hole for hole in holes])


    temp2 = delete_matrix_cols(temp1, del_col)

    println("stabilizer_holes: $(temp2)")

    return temp2
end

function get_stabilizer(R::Integer, C::Integer, X_nbr, Z_nbr)
    num_qubits = R*C

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

    # r: 1, 2, ..., d ;; c = 2
    sc_lx = [ 1+ r*C + 1 for r in 0:(R-1) ]

    for idx in sc_lx
        stabilizer[num_qubits, idx] = true
    end

    nq = num_qubits

    P = GF2.(stabilizer)
    rP = LinearAlgebra.rank(P)

    # println("\nA=$(A)")

    Q = rref(P)
    rQ = LinearAlgebra.rank(Q)

    println("rPQ: $(rP), $(rQ)")

    s_temp = stabilizer[1:num_qubits-1, :]

    A = GF2.(s_temp)
    rA = LinearAlgebra.rank(A)

    # println("\nA=$(A)")

    B = rref(A)
    rB = LinearAlgebra.rank(B)

    println("rAB: $(rA), $(rB)")

    a_temp = s_temp

    # for c in [13, 33, nq+13, nq+33]
    #     for r in 1:R
    #         a_temp[r,c] = false
    #     end
    # end

    # println("\nstabilizer=$(stabilizer)")


    # println("\na_temp=$(a_temp)")

    for q1 in 1:nq
        for q2 in 1:nq

            println("q1, q2 = $(q1), $(q2)")

            a_new = delete_matrix_cols(a_temp, [q1, q2, nq+q1, nq+q2]) 
            
            l_op = stabilizer[num_qubits, :]
           
            s_holes = get_logical_op_holes(R, C, [q1, q2], l_op, "Z")

            temp = s_holes
            a_new = vcat(a_new, temp)


            C = GF2.(a_new)
            D = rref(C)
            
            rD = LinearAlgebra.rank(D)
            println("q1=$(q1), q2=$(q2) Ranks-QD: $(rQ), $(rD)")


            if rD<rQ
                println("YAY: q1=$(q1), q2=$(q2) Ranks-QD: $(rQ), $(rD)")
            end
        end
    end

    return stabilizer

end

function check_surface_code_decoder(R::Integer, C::Integer)
    @info "Initialization Stage"
    t0 = time()
    begin
        (X_nbr, Z_nbr) = get_nbr(R, C)
        stabilizer = get_stabilizer(R, C, X_nbr, Z_nbr)

        phases = get_phases(R, C)
    end

    # println("Encoded stabilizer : $(stabilizer)")

    println("X idxs: $(X_nbr)")
    println("Z idxs: $(Z_nbr)")

    @info "Configuration Init"
    t1 = time()
    begin
        (ρ01, ϕ_x1, cfg1) = get_config(stabilizer, phases, X_nbr, Z_nbr, R, C)
        cfgs1 = QuantSymEx(cfg1)
    end

    @info "SMT constraint Generation"
    t2 = time()

    begin
        res = false
        for cfg in cfgs1
            if res && !generate_constraints(
                cfg.ρ, ρ01, (ϕ_x1 #=& ϕ_z2=#, cfg.ϕ[1], cfg.ϕ[2]),
                )
                res = false
                break
            end
        end

    end

    println("res after SymEx: $(res)")

    @info "SMT Solver Stage"
    t3 = time()

    begin
        if res
            res = solve_constraints(`bitwuzla --smt-comp-mode true -rwl 0 -S kissat`)
        end
    end

    t4 = time()

    res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
end

# check_surface_code_decoder(3) # precompile time
@info "precompile done..."

open("surface_code.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 3:3
        res_d, all, init, config, cons_gen, cons_sol = check_surface_code_decoder(9, 5)
        println("d,res,nq,all,init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(d*d),$(all),$(init),$(config),$(cons_sol),$(cons_gen)")
        println(io, "$(d),$(res_d),$(45),$(all),$(init),$(config),$(cons_sol),$(cons_gen)")
    end
end
