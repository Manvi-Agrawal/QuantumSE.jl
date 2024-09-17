module QEC_Helper

using QuantumSE

function get_stabilizer_from_nbr(num_qubits::Integer, X_nbr, Z_nbr, l_op)

    stabilizer = fill(false, num_qubits, 2*num_qubits)

    for row in 1:length(X_nbr)
        for idx in X_nbr[row]
            # println("Set $(row), $(idx) to true...")
            stabilizer[row, idx] = true
        end
    end
    
    for row in 1 : length(Z_nbr)
        for idx in Z_nbr[row]
            s_row = row + length(X_nbr)
            s_col = num_qubits+idx
            # println("Set $(s_row), $(s_col) to true...")
            stabilizer[s_row, s_col] = true
        end
    end 

    for row in 1:length(l_op)
        s_row = row + length(X_nbr) + length(Z_nbr)
        (l_op_row, s_type) = l_op[row]
        
        for idx in l_op_row
            s_col = s_type == "X" ? idx : num_qubits+idx
            stabilizer[num_qubits, s_col] = true
        end
    end

    return stabilizer

end




# export get_stabilizer

end