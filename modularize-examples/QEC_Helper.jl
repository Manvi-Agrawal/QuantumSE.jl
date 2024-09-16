module QEC_Helper

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

# export get_stabilizer

end