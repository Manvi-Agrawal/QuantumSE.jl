using QuantumSE
using Z3

using SimpleGF2
using LinearAlgebra
using InvertedIndices

function sanitise_nbr(nbr, holes)
    return map(x -> filter(e -> e∉holes, x), nbr) 
end

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

    # println("X 2q: $(x_2q)")


    z_2q_right = [ 1 + r*C+ (C-1) for r in 0:2:(R-2) ]
    z_2q_left = [ 1 + r*C for r in 1:2:(R-2) ]
    z_2q = vcat(z_2q_left, z_2q_right)

    #  println("Z 2q: $(z_2q)")

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

    return stabilizer

end

function get_holes(stabilizer, R::Integer, C::Integer, X_nbr, Z_nbr)

    P = GF2.(stabilizer)
    rP = LinearAlgebra.rank(P)

    Q = rref(P)
    rQ = LinearAlgebra.rank(Q)

    nq = R*C

    for hole in 1:nq
        a_temp = stabilizer[1:nq-1, :]
        a_new = delete_matrix_cols(a_temp, [hole, nq+hole]) 
        
        


        C = GF2.(a_new)
        rC = LinearAlgebra.rank(C)

        D = rref(C)
        
        rD = LinearAlgebra.rank(D)
        println("hole = $(hole), Ranks-CD: $(rC), $(rD)")


        if rD<rQ
            println("YAY: hole = $(hole), Ranks-QD: $(rQ), $(rD)")
        end

        # l_op = stabilizer[nq, :]
        
        # s_holes = get_logical_op_holes(R, C, [hole], l_op, "Z")

        # temp = s_holes
        # a_new = vcat(a_new, temp)
    end

end

function check_surface_code_decoder(R::Integer, C::Integer)
    @info "Initialization Stage"
    res = true

    println("d,q_h,res,nq,all,init,config,cons_gen,cons_sol")

    ts = time()

    begin
        num_qubits = R*C
        (X_nbr, Z_nbr) = get_nbr(R, C)
        # q_h = (R÷2)*C + C

        # println("q_h=$(q_h)")

        println("X_nbr: $(X_nbr)")
        println("Z_nbr: $(Z_nbr)")


        stabilizer = get_stabilizer(R, C, X_nbr, Z_nbr)
        # println("Encoded stabilizer : $(stabilizer)")

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

        println("rAB: $(rP), $(rQ)")

        get_holes(stabilizer, R, C, X_nbr, Z_nbr)
    end


    te = time()

    println("Total time: $(te-ts)")

    # res, t4-t0, t1-t0, t2-t1, t3-t2, t4-t3
    res, te-ts
end



check_surface_code_decoder(9, 5)
