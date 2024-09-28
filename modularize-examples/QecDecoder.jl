include("QEC_Decoder_Defaults.jl")
using .QEC_Decoder_Defaults

Base.@kwdef mutable struct QecDecoder
    d::Integer
    num_qubits::Integer
    ctx
    _xadj = missing
    _zadj = missing
    nx = missing
    nz = missing
    decoder = QEC_Decoder_Defaults.qec_decoder_ckt
    decoder_params = (nx, nz, num_qubits, d, ctx)
    x_syndrome_circuit = QEC_Decoder_Defaults.x_syndrome_circuit
    z_syndrome_circuit = QEC_Decoder_Defaults.z_syndrome_circuit
    decoder_algo_xz = QEC_Decoder_Defaults.decoder_algo_xz
    bug = QEC_Decoder_Defaults.bug
end