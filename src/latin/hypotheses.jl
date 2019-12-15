

function test2(fam :: Family)
    tmp = []
    orbFam = orb(fam)
    res = zeros(length(orbFam))
    for (i, ff) in enumerate(orbFam)
        total_assoc = assoc(ff)
        res[i] = sum(total_assoc) / length(total_assoc)
        push!(tmp, countmap(total_assoc))
    end 
    return res, tmp
end

function test3(fams :: Vector{Family{T}}) where T
    min_assoc = zeros(Int, length(fams))
    mean_assoc = zeros(length(fams))
    for (i, fam) in enumerate(fams)
        total_assoc = assoc(fam)
        min_assoc[i] = minimum(total_assoc)
        mean_assoc[i] = sum(total_assoc) / length(total_assoc)
    end
    return min_assoc, mean_assoc
end   


function test5(assoc4)
    res = [[] for i in 1 : 56]
    for i in 1 : 56
        res[i] = min_pairings(assoc4[1][i])
    end
    return res
end

function test6(assoc3)
    res = [[] for i in 1 : 2]
    for i in 1 : 2
        res[i] = min_pairings(assoc3[1][i])
    end
    return res
end


n = length(prs3)
for i in 1 : n
    prs = prs3[i]
    for pr in prs
        print(switch(pr) in prs)
        print()
    end
    println()
end

Family(
[[2,3], [2,4], [2,5], [3,4], [3,5], [4,5]],
[[1], [1,3], [1,4], [1,5], [3,4], [3,5], [4,5]],
[[1], [2], [1,2], [1,4], [1,5], [2,4], [2,5], [4,5]],
[[1], [2], [3], [1,2], [1,3], [1,5], [2,3], [2,5], [3,5]],
[[1], [2], [3], [4], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
)


Family(
[[2,3], [2,4], [3,4]],
[[1,3], [1,4], [3,4]],
[[1,2], [1,4], [2,4]],
[[1,2], [1,3], [2,3]]
)


function experimental_family(n)
    fam = []
    for fn in 1 : n
        fun = Vector{Vector{Int}}([])
        for i in 1 : n
            for j in 1 : (i - 1)
                if (i != fn) & (j != fn)
                    push!(fun, [i, j])
                end
            end
            if (i < fn)
                push!(fun, [i])
            end
        end
        push!(fam, fun)
    end
    return Family(fam...)
end
