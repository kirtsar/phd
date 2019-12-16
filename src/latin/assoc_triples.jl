# construct all latin squares of given size
# find the number of associative triples
import Base.getindex
import Base.show


struct AssocAnalyzer{T}
    pairings :: Vector{T}
    num_assoc :: Vector{Int}
end


function show(io :: IO, aa :: AssocAnalyzer)
    println("AssocAnalyzer of size $(length(aa.num_assoc))")
end


function assoc_init(n :: Int)
    pairs = generate_all_pairings(n)
    num_assoc = zeros(Int, length(pairs))
    return AssocAnalyzer(pairs, num_assoc)
end


function assoc_init(f :: Family)
    m = length(fam)
    return assoc_init(m)
end


function (a :: AssocAnalyzer)(f :: Family)
    n = 2^length(fam)
    table = zeros(Int, (n, n))
    ls = LatinSquare(table)
    for (i, p) in enumerate(a.pairings)
        a.num_assoc[i] = assoc!(ls, fam, p)
    end
end


function getindex(a :: AssocAnalyzer, i :: Int)
    return (a.pairings[i], a.num_assoc[i])
end


function getindex(a :: AssocAnalyzer, I) 
    return [a[i] for i in I]
end


function minimal(a :: AssocAnalyzer)
    i = argmin(a.num_assoc)
    return a[i]
end


function minimal(a :: AssocAnalyzer, i :: Int)
    inds = sortperm(a.num_assoc)[1 : i]
    return a[inds]
end



# find number of associative triples
# in latin square
function assoc(ls :: LatinSquare)
    m = length(ls)
    assocNum = 0
    for i in 0 : m - 1
        for j in 0 : m - 1
            for k in 0 : m - 1
                # (ij)k == i(jk) ?
                if ls(ls(i, j), k) == ls(i, ls(j, k))
                    assocNum += 1
                end
            end
        end
    end
    return assocNum
end


# return associative triples
function triples(ls :: LatinSquare)
    m = length(ls)
    assoc_triples = Vector{Tuple{Int, Int, Int}}([])
    for i in 0 : m - 1
        for j in 0 : m - 1
            for k in 0 : m - 1
                # (ij)k == i(jk) ?
                if ls(ls(i, j), k) == ls(i, ls(j, k))
                    #assocNum += 1
                    push!(assoc_triples, (i, j, k))
                end
            end
        end
    end
    return assoc_triples
end


# return triples for family and pairings
function triples(fam :: Family{T}, pi :: Pairing{T}) where T
    m = length(fam)
    ls = LatinSquare(fam, pi)
    assoc_triples = triples(ls)
    return assoc_triples
end



function assoc(fam :: Family)
    aa = assoc_init(fam)
    aa(fam)
    return aa
end


function assoc!(
        ls :: LatinSquare, 
        fam :: Family, 
        p :: Pairing)
    latin!(ls, fam, p)
    total_assoc = assoc(ls)
    return total_assoc
end

#=
function total_assoc(fams :: Vector{Family{T}}) where T
    min_assoc = zeros(Int, length(fams))
    mean_assoc = zeros(length(fams))
    m = length(fams[1])
    pairs = generate_all_pairings(m)
    n = 2^m - 1
    table = zeros(Int, (n + 1, n + 1))
    ls = LatinSquare(table)
    #fout = open("totalmin.txt", "w")
    total_assoc = zeros(Int, length(pairs))
    for (i, fam) in enumerate(fams)
        assoc!(ls, fam, pairs, total_assoc)
        min_assoc[i] = minimum(total_assoc)
        mean_assoc[i] = sum(total_assoc) / length(total_assoc)
        #print(fout, minimum(total_assoc))
        #print(fout, ",")
        #println(fout, sum(total_assoc) / length(total_assoc))
    end
    #close(fout)
    return min_assoc, mean_assoc
end  


function save_assoc(trans, tmin, tmean; filename = "assoc4.clf")
    n = length(trans)
    familyLen = length(trans[1])
    out = open("DATA/" * filename, "w")
    println(out, familyLen)
    println(out, n)
    for i in 1 : n
        #### PRINT WHOLE FAMILY ####
        curr_family = trans[i]
        # print current family index
        println(out, i)
        zhegArray = curr_family.fs
        for j in 1 : familyLen
            currZheg = zhegArray[j]
            currANF = currZheg.ANF
            for k in 1 : length(currANF)
                currMonom = currANF[k]
                val = currMonom.val
                # print each monom in file
                print(out, val)
                print(out, " ")
            end
            # between different families
            println(out, "")
        end
        #### PRINT tmin #####
        println(out, tmin[i])
        #### PRINT tmean ####
        println(out, tmean[i])
    end
    close(out)
end


function load_assoc(filename)
    in = open("../data/" * filename, "r")
    familyLen = parse(Int, readline(in))
    transLen = parse(Int, readline(in))
    tmin = zeros(Int, transLen)
    tmean = zeros(transLen)
    trans = Vector{Family{Int}}([])
    for j in 1 : transLen
        familyNumber = parse(Int, readline(in))
        @assert j == familyNumber
        zhegArray = []
        for k in 1 : familyLen
            monoms = parse.(Int, split(readline(in)))
            zhegfun = ZhegFun(monoms)
            push!(zhegArray, zhegfun)
        end
        fam = Family(zhegArray...)
        push!(trans, fam)
        tmin[j] = parse(Int, readline(in))
        tmean[j] = parse(Float64, readline(in))
    end
    close(in)
    return trans, tmin, tmean
end

=#

#=
function resave_assoc(filename)
    assocClf = load_assoc(filename)
    assocClf = sort_assoc(assocClf)
    m = length(assocClf[1])
    for i in 1 : m
        assocClf[1][i] = minimal_repr(assocClf[1][i])
    end
    tmp = sortslices([assocClf[1] assocClf[2] assocClf[3]], 
                      by = x -> (x[2], x[3]),
                      dims = 1)
    trans = tmp[:, 1]
    tmin = tmp[:, 2]
    tmean = tmp[:, 3]
    save_assoc(trans, tmin, tmean; filename = filename)
end
=#