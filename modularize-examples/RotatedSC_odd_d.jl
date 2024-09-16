using QuantumSE

include("QEC_Defaults.jl")
using .QEC_Defaults

include("QEC_Pipeline.jl")
using .QEC_Pipeline

include("QEC_Helper.jl")
using .QEC_Helper

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

function decoder_bug(ρ, r_x, r_z, d)
    println("RSC_bug")
    e = reduce(&, r_z[1:(d-1)÷2])

    # sX(ρ, 1, e)
    sX(1, e)

    println("RSC_bug end")
    nothing

end


# rsc_d3 = QEC_Pipeline.QecDecoderConfig(
#             d=3,
#             X_nbr=get_nbr(3)[1],
#             Z_nbr=get_nbr(3)[2],
#             phases=get_phases(3),
#             ctx=ctx)

# QEC_Pipeline.check_qec_decoder(rsc_d3) # precompile time
# @info "precompile done..."



open("rsc.csv", "w") do io
    println(io, "d,res,nq,all,init,config,cons_gen,cons_sol")

    for d in 3:2:7
        tm2 = time()
        nq = d*d
        
        (X_nbr, Z_nbr) = get_nbr(d)

        stabilizer = QEC_Helper.get_stabilizer(d, X_nbr, Z_nbr)

        _xadj(j) = X_nbr[j]
        _zadj(j) = Z_nbr[j]

        rsc_decoder = QEC_Helper.qec_decoder
        xlim = (nq-1)÷2
        zlim = (nq-1)÷2
        rsc_bug = QEC_Defaults.bug
        rsc_decoder_params = (ctx, d, nq, xlim, zlim, _xadj, _zadj, rsc_bug)

        
        rsc_decoder_config = QEC_Pipeline.QecDecoderConfig(
            d=d,
            _xadj=_xadj,
            _zadj=_zadj,
            stabilizer=stabilizer,
            phases= get_phases(d),
            ctx=ctx,
            bug = rsc_bug,
            decoder=rsc_decoder,
            decoder_params=rsc_decoder_params)

        tm1 = time()
        

        # println("X_nbr: $(X_nbr)")
        # println("Z_nbr: $(Z_nbr)")

        # println("Phases: $(phases)")

        

        res_d, all, init, config, cons_gen, cons_sol = QEC_Pipeline.check_qec_decoder(rsc_decoder_config)

        init_config = (tm2-tm1)
        all += init_config

        println("d,res,nq,all,init_config, init,config,cons_gen,cons_sol")
        println("$(d),$(res_d),$(d*d),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
        println(io, "$(d),$(res_d),$(d*d),$(all),$(init_config),$(init),$(config),$(cons_gen),$(cons_sol)")
    end
end
