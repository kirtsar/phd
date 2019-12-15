import Base.length
import Base.getindex

# family of zhegalkin functions 

@auto_hash_equals mutable struct Family{T}
    fs :: Vector{ZhegFun{T}}
end


function Family(arr :: Vector{Vector{T}}) where T
    return Family(ZhegFun.(arr))
end


function Family(fargs...)
    return Family([fargs...])
end


# apply family of functions
function (funs :: Family)(x)
    s = 0
    for i in length(funs.fs) : -1 : 1
        s <<= 1
        s |= funs.fs[i](x)
    end
    return s
end


# return size of the family

function length(fam :: Family)
    return length(fam.fs)
end


# get i'th function
function getindex(fam :: Family, i)
    return fam.fs[i]
end