module QEC_Helper

using QuantumSE

function encoder_from_nbr(num_qubits::Integer, X_nbr, Z_nbr, l_op)
    stabilizer = fill(false, num_qubits, 2*num_qubits)
    phases = Vector{Z3.ExprAllocated}(undef, num_qubits)

    lx = _bv_const(ctx, "lx")
    lz = _bv_const(ctx, "lz")

    for row in 1:length(X_nbr)
        for idx in X_nbr[row]
            # println("Set $(row), $(idx) to true...")
            stabilizer[row, idx] = true
            phases[row] = _bv_val(ctx, 0)
        end
    end
    
    for row in 1 : length(Z_nbr)
        for idx in Z_nbr[row]
            s_row = row + length(X_nbr)
            s_col = num_qubits+idx
            # println("Set $(s_row), $(s_col) to true...")
            stabilizer[s_row, s_col] = true
            phases[s_row] = _bv_val(ctx, 0)
        end
    end 

    for row in 1:length(l_op)
        s_row = row + length(X_nbr) + length(Z_nbr)
        (l_op_row, s_type) = l_op[row]
        
        for idx in l_op_row
            s_row = row + length(X_nbr)+length(Z_nbr)
            s_col = s_type == "X" ? idx : num_qubits+idx
            stabilizer[s_row, s_col] = true
            phases[s_row] = s_type == "X" ? lx : lz
        end
    end

    return (stabilizer, phases)
end




# export get_stabilizer

end