function _zadj(d, q_idx)
    if q_idx == 1
        return [3]
    elseif q_idx == 3
        return [1, 4]
    elseif q_idx == 4
        return [3]
    else
        return []
    end
end

res= _zadj(2, 1)
println("Res: $(res[1])")