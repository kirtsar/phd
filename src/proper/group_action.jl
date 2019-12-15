using Combinatorics

function to_bits!(m :: Monom, bits)
    curr = m.val
    i = 1
    for i in 1 : length(bits)
        bits[i] = curr & 0x1
        curr >>= 1
    end
    return nothing
end

function from_b(bits)
    res = 0
    for i in length(bits) : -1 : 1
        res <<= 1
        res |= bits[i]
    end
    return Monom(res)
end


# permute monom
function permute(m :: Monom, perm)
    n = length(perm)
    bits = zeros(UInt8, n)
    to_bits!(m, bits)
    return from_b(bits[perm])
end



# permute proper family of functions
# fi(x1, ... xn) -> f_pi(i)(x_pi(1), ... x_pi(n))
# this family is also proper
# matrix graph is the same
function permute(f :: ZhegFun, perm)
    g = ZhegFun(Monom.(zeros(Int, length(f.ANF))))
    for (ind, mon) in enumerate(f.ANF)
        g.ANF[ind] = permute(mon, perm)
    end
    sort!(g.ANF, by = x -> x.val)
    return g
end


# given a family of functions
# apply all possible permutations to this family
# fi(x1, ... xn) -> f_pi(i)(x_pi(1), ... x_pi(n))
# obtain orbit of the family
function orb(fam :: Family{T}) where T
    n = length(fam.fs)
    perms = permutations(collect(1 : n))
    resOrbit = Set{Family{T}}([])
    for perm in perms
        funArray = Vector{ZhegFun{T}}([])
        for i in 1 : n
            fi = permute(fam.fs[i], perm)
            push!(funArray, fi)
        end
        push!(resOrbit, Family(funArray[perm]...))
    end
    return resOrbit
end


# obtain stabilizer of the family
# under action fi(x1, ... xn) -> f_pi(i)(x_pi(1), ... x_pi(n))

function stab(fam :: Family{T}) where T
    n = length(fam)
    perms = permutations(collect(1 : n))
    stabilizers = Vector{Vector{Int}}([])
    for perm in perms
        funArray = Vector{ZhegFun{T}}([])
        for i in 1 : n
            fi = permute(fam.fs[i], perm)
            push!(funArray, fi)
        end
        fam_new = Family(funArray[perm]...)
        if fam_new == fam 
            push!(stabilizers, perm)
        end
    end
    return stabilizers
end


# represent number as vector of 0-1
function as_vec!(funVec, numb :: T) where T <: Integer
    n = length(funVec)
    for i in n : -1 : 1
        funVec[i] = numb & 0x1
        numb >>= 1
    end
    return nothing
end


# obtain orbit under the action 
# fi(x1, ... xn) -> fi(x1, ... xn) + const_i
function additive_orb(fam :: Family{T}) where T
    n = length(fam)
    shiftOrb = Vector{Family{T}}()
    currFam = deepcopy(fam)
    shiftVec = zeros(Int, n)
    for shift in 0 : 2^n - 1
        as_vec!(shiftVec, shift)
        for j in 1 : n
            if (shiftVec[j] == 1) && !(Monom(0) in currFam[j].ANF)
                push!(currFam[j].ANF, Monom(0))
                sort!(currFam[j].ANF, by = (x -> x.val))
            elseif (shiftVec[j] == 0) && (Monom(0) in currFam[j].ANF)
                filter!(x -> x != Monom(0), currFam[j].ANF)
                sort!(currFam[j].ANF, by = (x -> x.val))
            end
        end
        push!(shiftOrb, deepcopy(currFam))
    end
    return shiftOrb
end