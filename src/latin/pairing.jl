@auto_hash_equals mutable struct Pairing{T}
    funs :: Family{T}
end


function length(pr :: Pairing{T}) where T
    return length(pr.funs)
end


function get_bit_pair(x, y)
    res = ((y & 0x1) << 1) + (x & 0x1)
    return res 
end


function (pair :: Pairing)(x, y)
    res = 0
    base = 1
    n = length(pair.funs)
    for i in 1 : n
        # bitpair is (x_i, y_i)
        bitpair = get_bit_pair(x, y)
        current_bit = pair.funs[i](bitpair)
        res += current_bit * base
        base <<= 1
        x >>= 1
        y >>= 1
    end
    return res
end

function generate_all_pairings(numVars)
    pairSet = Vector{Pairing}([])
    # only possible coefs for function f(x, y) is 0, 1, 2, 3
    famReprSet = Set()
    for i in 1 : numVars
        push!(famReprSet, powerset([0, 1, 2, 3]))
        #push!(famReprSet, powerset([1, 2, 3]))
    end

    for famRepr in product(famReprSet...)
        pair = Pairing(Family(famRepr...))
        push!(pairSet, pair)
    end

    return pairSet
end

# NB : pairing is slow compared to 
# direct usage of function : 
# pi(x, y) => 32.4 ns
# fst(x, y) => 0.023 ns
# with parametrization:
# pi(x, y) => 15 ns



# example test
function test_pairings()
    # [-1] or [] corresponds to 0
    # [0] corresponds to 1
    # [1] corresponds to x
    # [2] corresponds to y
    # [3] corresponds to xy
    # [1, 2] corresponds to (x XOR y)
    
    zeroZh = Int[]
    oneZh = [0]
    # we must obtain 6 in every arg
    pair = Pairing(Family(ZhegFun.([zeroZh, oneZh, oneZh])))
    test1 = true
    for i in 0 : 7
        for j in 0 : 7 
            if pair(i, j) != 6
                test1 = false
            end
        end
    end

    if test1
        println("test1 passed")
    end

    # now we must obtain first argument
    pair = Pairing(Family(ZhegFun.([[1], [1], [1]])))
    test2 = true
    for i in 0 : 7
        for j in 0 : 7 
            if (pair(i, j) != i)
                test2 = false
            end
        end
    end
    if test2
        println("test2 passed")
    end

    # here we must obtain i XOR j
    pair = Pairing(Family(ZhegFun.([[1, 2], [1, 2], [1, 2]])))
    test3 = true
    for i in 0 : 7
        for j in 0 : 7 
            if (pair(i, j) != xor(i, j))
                test3 = false
            end
        end
    end

    if test3
        println("test3 passed")
    end
end


function permute(pr :: Pairing{T}, perm) where T
    return Pairing(Family(pr.funs[perm]))
end


function orb(pr :: Pairing{T}) where T
    n = length(pr)
    res = Vector{Pairing{T}}([])
    for perm in permutations(1:n)
        prnext = permute(pr, perm)
        push!(res, prnext)
    end
    return res
end


function min_pairings(fam :: Family{T}) where T
    m = length(fam)
    pairs = generate_all_pairings(m)
    assoc_number = assoc(fam)
    min_assoc = minimum(assoc_number)
    return pairs[assoc_number .== min_assoc]
end

function switch(pr :: Pairing{T}) where T
    newfuns = Vector{ZhegFun{Int}}([])
    n = length(pr.funs)
    for i in 1 : n
        fun = pr.funs[i]
        push!(newfuns, permute(fun, [2, 1]))
    end
    return Pairing(Family(newfuns))
end