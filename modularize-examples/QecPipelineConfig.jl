include("QEC_Defaults.jl")
using .QEC_Defaults

Base.@kwdef mutable struct QecPipelineConfig
    # d::Integer
    # num_qubits::Integer
    stabilizer
    phases
    # ctx
    # _xadj = missing
    # _zadj = missing
    # nx = missing
    # nz = missing
    decoder = QEC_Defaults.qec_decoder()
end
