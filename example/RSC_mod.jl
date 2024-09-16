using QuantumSE


include("QEC_Pipeline.jl")
using .QEC_Pipeline

using Z3

N4 = 4
N2 = 2

AD4 = 4
AD2_H = 2
AD2_V = 2

ctx = Context()


function get_nbr(d::Integer)
    X_idxs = [[2, 3, 5, 6], [4, 5, 7, 8], [1, 2], [8, 9]]

    Z_idxs = [[1, 2, 4, 5], [5, 6, 8, 9], [4, 7], [3, 6]]

    return (X_idxs, Z_idxs)

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



# check_qec_decoder(3) # precompile time
# @info "precompile done..."

open("surface_code.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 3:2:3
        (X_nbr, Z_nbr) = get_nbr(d)
        phases = get_phases(d)
        println("X_nbr: $(X_nbr)")
        println("Z_nbr: $(Z_nbr)")

        println("Phases: $(phases)")

        res_d, all, init, config, cons_gen, cons_sol = QEC_Pipeline.check_qec_decoder(d, X_nbr, Z_nbr, phases, ctx)
        println("d,res,nq,all,init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
        println(io, "$(d),$(res_d),$(d*d),$(all),$(init),$(config),$(cons_gen),$(cons_sol)")
    end
end
