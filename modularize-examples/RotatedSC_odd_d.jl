using QuantumSE

include("QEC_Decoder_Defaults.jl")
using .QEC_Decoder_Defaults

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

function rsc_l_op(d::Integer)
    # r: 1, 2, ..., d ;; c = 1
    sc_lx = [ 1+ r*d + 1 for r in 0:(d-1) ]
    return [(sc_lx, "X")]
end

function get_rsc_decoder_config(d::Integer)
    nq = d*d
        
    (X_nbr, Z_nbr) = get_nbr(d)
    (stabilizer,phases) = QEC_Helper.encoder_from_nbr(nq, X_nbr, Z_nbr, rsc_l_op(d))
    _xadj(j) = X_nbr[j]
    _zadj(j) = Z_nbr[j]

    rsc_decoder_config = QEC_Pipeline.QecPipelineConfig(
        d=d,
        num_qubits = nq,
        _xadj=_xadj,
        _zadj=_zadj,
        nx=(nq-1)÷2,
        nz=(nq-1)÷2,
        stabilizer=stabilizer,
        phases= get_phases(nq),
        ctx=ctx)

    return rsc_decoder_config
end

QEC_Pipeline.qec_runner("rsc.csv",get_rsc_decoder_config, 3:2:11)
