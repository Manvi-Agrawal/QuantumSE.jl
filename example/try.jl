# function modify_stabilizer(stabilizer)
#     stabilizer[1, 1] = true
#     println("inside modify_stabilizer: $(stabilizer)")
# end

# stabilizer = fill(false, 4, 5)
# modify_stabilizer(stabilizer)
# println("outside modify_stabilizer: $(stabilizer)") # Verdict-nodified

### Ranges
# a=1:2:9
# [t for t in a]

# function matrix_to_vector_of_vectors(matrix)
#     vv = []
#     for row in 1:3
#         push!(vv, matrix[row, :])
#     end
#     return vv
# end

# # Example usage:
# matrix = [1 2 3; 4 5 6; 7 8 9]
# vector_of_vectors = matrix_to_vector_of_vectors(matrix)
# println(vector_of_vectors)


using LinearAlgebra

using SimpleGF2

# function remove_zero_rows(a)
#     # Filter rows that have at least one non-zero element
#     return a[vec(mapslices(col -> any(col .!= 0), a, dims = 2)), :]
# end

function remove_false_rows(a)
    # Filter rows that have at least one non-zero element
    return a[vec(mapslices(col -> any(col .!= false), a, dims = 2)), :]
end

# X_nbr=[[2, 3, 7, 8], [4, 5, 9, 10],  [42, 43], [44, 45]]
# Z_nbr=[[1, 2, 6, 7], [3, 4, 8, 9], [7, 8, 12, 13], [9, 10, 14, 15], [11, 12, 16, 17], [13, 14, 18, 19], [17, 18, 22, 23], [19, 20, 24], [21, 22, 26, 27], [23, 24, 28, 29], [27, 28, 32, 33], [29, 30, 34, 35], [31, 32, 36, 37], [33, 34, 38, 39], [37, 38, 42, 43], [39, 40, 44, 45], [6, 15], [16], [26, 35], [36, 45], [5, 14], [15, 24], [34], [35, 44]]

# A = [1 1 0 0 0 0 0 0; 0 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0]
A = [1 1 0 0 ; 0 1 1 0 ; 0 0 0 0]

B = GF2.(A)
C = [Int64(x.val)==1 for x in B]

C_n = remove_false_rows(C)
println(C_n)

function populate_neighbour(stabilizer)
    nR = length(stabilizer[:, 1])
    nC = length(stabilizer[1, :]) ÷ 2

    println("[nR, nC]: [$(nR), $(nC)]")

    X_nbr = []

    for r in 1:nR
        nl = []
        for c in 1:nC
            if stabilizer[r, c]
                val = c
                println("Pushing $(val):: at $(r), $(c)")

                push!(nl, val)
            end
        end
        push!(X_nbr, nl)
    end

    return X_nbr

end
# row = [ 0 0 1 1]
# a_new = vcat(matrix, row)


