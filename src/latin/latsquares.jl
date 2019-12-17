import Base.length


struct LatinSquare{Int}
    sq :: Matrix{Int}
end


function (ls :: LatinSquare)(i :: Int, j :: Int)
    return ls.sq[i + 1, j + 1]
end


function length(ls :: LatinSquare)
    return size(ls.sq)[1]
end

# construct latin square from given 
# family and pairing
function latin(fam :: Family, pair :: Pairing)
    n = length(fam)
    m = 2^n - 1
    table = zeros(Int, (m + 1, m + 1))
    ls = LatinSquare(table)
    latin!(ls, fam, pair)
    return LatinSquare(table)
end


# in-place constructor of latin square
function latin!(LS :: LatinSquare, fam :: Family, pair :: Pairing)
    n = length(fam)
    m = 2^n - 1
    # table = zeros(Int, (m + 1, m + 1))
    for i in 0 : m
        for j in 0 : m
            # (i, j) -> i + j + f o pi (i, j)
            LS.sq[i + 1, j + 1] = xor(xor(i, j), fam(pair(i, j)))
        end
    end
    # return LatinSquare(table)
end


function test_LS()
    zeroZh = Int[]
    fam = Family(ZhegFun.([zeroZh, zeroZh, zeroZh]))
    # does not depend on pairing
    pair = Pairing(Family(ZhegFun.([[1, 2], [1, 2], [1, 2]])))
    LS = latin(fam, pair)
    return LS.sq
end



    