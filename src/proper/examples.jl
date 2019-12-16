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


# generates ortho family of size n
function ortho_family(n)
    mon = vcat(collect(2 : n), collect(1 : n))
    funs = Vector{ZhegFun}([])
    for i in 1 : n
        m1 = monom(mon[i : i+(n-2)])
        m2 = monom(mon[i+1 : i+(n-2)])
        push!(funs, ZhegFun([m1, m2]))
    end
    return Family(funs...)
end