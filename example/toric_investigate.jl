using QuantumSE
using Z3

using SimpleGF2
using LinearAlgebra
using InvertedIndices

ctx = Context()

function _xadj(d, idx)
    ai = (idx-1)÷d
    bi = (idx-1)%d

    [d*d+ai*d+bi, ai*d+(bi+d-1)%d, ai*d+bi, d*d+((ai+d-1)%d)*d+bi] .+ 1
end

function _zadj(d, idx)
    ai = (idx-1)÷d
    bi = (idx-1)%d

    [((ai+1)%d)*d+bi, d*d+ai*d+bi, d*d+ai*d+(bi+1)%d, ai*d+bi] .+ 1
end



function toric_x_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in _xadj(d, idx)
		s[j] = 1
	end

	s
end



function toric_z_s(d::Integer, idx::Integer)
	s = zeros(Bool, 4*d*d)

	for j in (_zadj(d, idx) .+ 2*d*d)
		s[j] = 1
	end

	s
end


function toric_lx1(d::Integer)
	s = zeros(Bool, 4*d*d)

	for i in 1:d
		s[i*d] = 1
	end

	s
end

function toric_lz1(d::Integer)
	s = zeros(Bool, 4*d*d)

	for i in 1:d
		s[3*d*d+i*d] = 1
	end

	s
end

function get_stabilizer(d)
    num_qubits = 2*d*d

    stabilizer = fill(0, num_qubits, 2*num_qubits)

    for i in 1:d*d-1
        stabilizer[i,:] = toric_x_s(d, i)
        stabilizer[i+d*d,:] = toric_z_s(d, i)
    end

    stabilizer[d*d,:] = toric_lx1(d)
    stabilizer[2*d*d,:] = toric_lz1(d)

    return stabilizer
end

function get_syndrome_stabilizer(d)
    num_qubits = 2*d*d

    stabilizer = fill(0, num_qubits-2, 2*num_qubits)

    for i in 1:d*d-1
        stabilizer[i,:] = toric_x_s(d, i)
        stabilizer[i+d*d-1,:] = toric_z_s(d, i)
    end

    return stabilizer
end

function check_toric_decoder(d::Integer)
    s = get_stabilizer(d)
    ss = get_syndrome_stabilizer(d)

    A = GF2.(s)
    rA = LinearAlgebra.rank(A)
    B = rref(A)
    rB = LinearAlgebra.rank(B)

    println("d=$(d), rA: $(rA), rB: $(rB)")
    # println("d=$(d), rA: $(rA)")



    C = GF2.(ss)
    rC = LinearAlgebra.rank(C)
    D = rref(C)
    rD = LinearAlgebra.rank(D)

    println("d=$(d), rC: $(rC), rD: $(rD)")
    # println("d=$(d), rC: $(rC)")


end





for d in 3:2:29
    check_toric_decoder(d)
end
