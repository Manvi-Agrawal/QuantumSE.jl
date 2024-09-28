module QEC_Defaults

function decoder()
    return QecDecoder(
        decoder = QEC_Decoder_Defaults.qec_decoder_ckt,
        decoder_params = (),
        x_syndrome_circuit = QEC_Decoder_Defaults.x_syndrome_circuit,
        z_syndrome_circuit = QEC_Decoder_Defaults.z_syndrome_circuit,
        decoder_algo_xz = QEC_Decoder_Defaults.decoder_algo_xz,
        bug = QEC_Decoder_Defaults.bug
    )

end

end