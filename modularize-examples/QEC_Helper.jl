module QEC_Helper

using QuantumSE

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

@qprog qec_decoder (nx, nz, nq, d, ctx) begin
    # println("Decoder start")

    s_x = [x_syndrome_circuit(j) for j in 1:nx]
    s_z = [z_syndrome_circuit(j) for j in 1:nz]

    r_x = decoder_algo_xz(ctx, d, s_x, "X", nq, _xadj)
    r_z = decoder_algo_xz(ctx, d, s_z, "Z", nq, _zadj)

    for j in 1:nq
        sZ(j, r_x[j])
        sX(j, r_z[j])
    end

    # a strange bug
    # bug(ρ, r_x, r_z, d)
    e = reduce(&, r_z[1:((d-1)÷2)])
    sX(1, e)

    # println("Decoder end")
end


# export get_stabilizer

end