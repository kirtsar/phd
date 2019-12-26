import Base.*
import Base.inv
import Base.show
import Base.getindex
import Base.setindex!

using OffsetArrays

struct NMap{T}
    img :: T
end


function NMap(T :: DataType, lens...)
    k = length(lens)
    x = zeros(T, lens...)
    img = OffsetArray(x, Tuple(-ones(Int, k)))
    return NMap(img)
end


function NMap(lens...)
    return NMap(Int, lens...)
end


function NMap(t :: Tuple)
    return NMap(Int, t...)
end


function ftype(nm :: NMap)
    return size(nm.img)
end


function Base.getindex(nm :: NMap, inds...)
    getindex(nm.img, inds...)
end


function Base.setindex!(nm :: NMap, val, inds...)
    setindex!(nm.img, val, inds...)
end


function Base.show(io :: IO, nm :: NMap)
    ft = ftype(nm)
    G = Grp(ft)
    maptype = repr(G)
    print(io, "Map $maptype -> Z")
end


function (nm :: NMap)(x...)
    return nm[x...]
end


function (nm :: NMap)(x :: Tuple)
    return nm(x...)
end


# represent map as table
# (i1, i2, ... , ik) -> f(i1, .., ik)
function view(nm :: NMap)
    ft = ftype(nm)
    G = Grp(ft)
    xs = elements(G)
    for x in xs
        println(x, " -> ", nm(x))
    end
end


struct NMapCollector{T}
    iter :: T
    ft :: Int
    st :: Int
    maxval :: Int
end


function NMapCollector(ft :: Tuple)
    iter = NMap(ft)
    st = zeros(Int, length(ft))
    maxval = prod(ft)
    return NMapCollector(iter, ft, st, maxval)
end

#=
function next!(nm :: NMapCollector)
    digits(nm.st, base = )
    nm.st += 1
end
=#




