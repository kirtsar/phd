BitType = Union{UInt8, UInt16, UInt32, Int}
using AutoHashEquals  # checking for equality

################################
# structures
# monoms in zhegalkin polys
struct Monom{T}
    val :: T
end


function monom(arr)
    res = 0
    for x in arr
        res += 2^(x-1)
    end
    return Monom(res)
end


# zhegalkin poly 
@auto_hash_equals mutable struct ZhegFun{T}
    ANF :: Vector{Monom{T}}
end


# zhegalkin poly from array
function ZhegFun(arr :: Vector{T}) where T <: BitType
    return ZhegFun(Monom.(arr))
end


function ZhegFun(monArr :: Vector{Vector{Int}})
    mons = zeros(Int, length(monArr))
    for (i, monRepr) in enumerate(monArr)
        val = 0
        for bit in monRepr
            val += 2^(bit - 1)
        end
        mons[i] = val
    end
    sort!(mons)
    return ZhegFun(mons)
end



################################
# APPLICATION OF DIFFERENT TYPES
# apply monom
function (m :: Monom)(x)
    return (m.val & x == m.val)
end


# apply zhegalkin poly

function (f :: ZhegFun)(x :: T) where T <: BitType
    res = 0
    for m in f.ANF
        res = xor(res, m(x))
    end
    return res
end



# find all essential variables of Zhegalkin Poly
function get_essentials(f :: ZhegFun)
    essentials = 0
    for coef in f.ANF
        essentials |= coef.val
    end
    return nzBits(essentials)
end






