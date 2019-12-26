maps1 = MFamilyCollector((2,3))
maps2 = MFamilyCollector((2,3))


# find all proper families on Z2 x Z3
function all_proper_23()
    propers = Any[]
    ft = (2,3)
    nm1 = NMap(2,3)
    nm2 = NMap(2,3)
    G = Grp(ft)
    inds1 = collect.(elements(G))
    inds2 = deepcopy(inds1)
    for fun1 in product(0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
        # initialize nm1 
        for (i, ind1) in enumerate(inds1)
            nm1[ind1...] = fun1[i]
        end
        for fun2 in product(0:2, 0:2, 0:2, 0:2, 0:2, 0:2)
            # initialize nm2
            for (j, ind2) in enumerate(inds2)
                nm2[ind2...] = fun2[j]
            end
            dfam = DFamily([nm1, nm2])
            if is_proper(dfam)
                push!(propers, deepcopy(dfam))
            end
        end
    end
    return propers
end


function all_proper_33()
    propers = Any[]
    ft = (3,3)
    nm1 = NMap(3,3)
    nm2 = NMap(3,3)
    G = Grp(ft)
    inds1 = collect.(elements(G))
    inds2 = deepcopy(inds1)
    for fun1 in tqdm(product(0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2))
        # initialize nm1 
        for (i, ind1) in enumerate(inds1)
            nm1[ind1...] = fun1[i]
        end
        for fun2 in product(0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2, 0:2)
            # initialize nm2
            for (j, ind2) in enumerate(inds2)
                nm2[ind2...] = fun2[j]
            end
            dfam = DFamily([nm1, nm2])
            if is_proper(dfam)
                push!(propers, deepcopy(dfam))
            end
        end
    end
    return propers
end


function all_proper_22()
    propers = Any[]
    ft = (2,2)
    nm1 = NMap(2,2)
    nm2 = NMap(2,2)
    G = Grp(ft)
    inds1 = collect.(elements(G))
    inds2 = deepcopy(inds1)
    for fun1 in product(0:1, 0:1, 0:1, 0:1)
        # initialize nm1 
        for (i, ind1) in enumerate(inds1)
            nm1[ind1...] = fun1[i]
        end
        for fun2 in product(0:1, 0:1, 0:1, 0:1)
            # initialize nm2
            for (j, ind2) in enumerate(inds2)
                nm2[ind2...] = fun2[j]
            end
            dfam = DFamily([nm1, nm2])
            if is_proper(dfam)
                push!(propers, deepcopy(dfam))
            end
        end
    end
    return propers
end



function all_proper_222()
    propers = Any[]
    ft = (2,2,2)
    nm1 = NMap(2,2,2)
    nm2 = NMap(2,2,2)
    nm3 = NMap(2,2,2)
    G = Grp(ft)
    inds1 = collect.(elements(G))
    inds2 = deepcopy(inds1)
    inds3 = deepcopy(inds1)
    for fun1 in product(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
        # initialize nm1 
        for (i, ind1) in enumerate(inds1)
            nm1[ind1...] = fun1[i]
        end
        for fun2 in product(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
            # initialize nm2
            for (j, ind2) in enumerate(inds2)
                nm2[ind2...] = fun2[j]
            end
            for fun3 in product(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
            # initialize nm2
                for (k, ind3) in enumerate(inds3)
                    nm3[ind3...] = fun3[k]
                end
                dfam = DFamily([nm1, nm2, nm3])
                if is_proper(dfam)
                    push!(propers, deepcopy(dfam))
                end
            end
        end
    end
    return propers
end




function reduced_33()
    propers = Any[]
    ft = (3,3)
    f1 = NMap(3,3)
    f2 = NMap(3,3)
    for fun1 in tqdm(Iterators.product(0:2, 0:2, 0:2))
        # initialize nm1 
        f1[0, 0] = f1[1, 0] = f1[2, 0] = fun1[1]
        f1[0, 1] = f1[1, 1] = f1[2, 1] = fun1[2]
        f1[0, 2] = f1[1, 2] = f1[2, 2] = fun1[3]
        for fun2 in Iterators.product(0:2, 0:2, 0:2)
            # initialize nm2
            f2[0, 0] = f2[0, 1] = f2[0, 2] = fun2[1]
            f2[1, 0] = f2[1, 1] = f2[1, 2] = fun2[2]
            f2[2, 0] = f2[2, 1] = f2[2, 2] = fun2[3]
            dfam = DFamily([f1, f2])
            if is_proper(dfam)
                push!(propers, deepcopy(dfam))
            end
        end
    end
    propers = Vector{typeof(propers[1])}(propers)
    return propers
end


# 2D-search ONLY!!!! 
function reduced_search(ft :: NTuple{2, T}) where T
    propers = []
    f1 = NMap(ft)
    f2 = NMap(ft)
    for fun1 in tqdm(Iterators.product(fill(0 : ft[1] - 1, ft[2])...))
        # initialize nm1 
        for j in 1 : ft[2]
            for i in 1 : ft[1]
                f1[i-1, j-1] = fun1[j]
            end
        end
        for fun2 in Iterators.product(fill(0 : ft[2] - 1, ft[1])...)
            # initialize nm2
            for i in 1 : ft[1]
                for j in 1 : ft[2]
                    f2[i-1, j-1] = fun2[i]
                end
            end
            dfam = DFamily([f1, f2])
            if is_proper(dfam)
                push!(propers, deepcopy(dfam))
            end
        end
    end
    propers = Vector{typeof(propers[1])}(propers)
    return propers
end
