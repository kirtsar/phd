import Base.length
import Base.getindex

using AutoHashEquals

# family of discrete functions 

@auto_hash_equals struct DFamily{T, S}
    fs :: Vector{T}
    ftype :: S
end


function DFamily(arr :: Vector{NMap{T}}) where T
    @assert all(ftype.(arr) .== [ftype(arr[1])])
    return DFamily(arr, ftype(arr[1]))
end


# return size of the family
function Base.length(fam :: DFamily)
    return length(ftype(fam))
end


# get i'th map
function Base.getindex(fam :: DFamily, i)
    return fam.fs[i]
end


# apply family of functions
function (fam :: DFamily)(x...)
    res = zeros(Int, length(fam))
    for i in 1 : length(fam)
        res[i] = fam[i](x...)
    end
    return Tuple(res)
end


function (fam :: DFamily)(x :: Tuple)
    return fam(x...)
end


function ftype(fam :: DFamily)
    return fam.ftype
end


function is_proper(fam :: DFamily)
    size = length(fam)
    ft = ftype(fam)
    k = length(ft)
    indices = CartesianIndices(ft) .- CartesianIndex(ones(Int, k)...)
    for i in indices
        for j in indices
            x = Tuple(i)
            y = Tuple(j)
            if x != y
                fx = fam(x)
                fy = fam(y)
                if !any((&).((x .!= y), (fx .== fy)))
                    return false
                end
            end
        end
    end
    return true
end



function Base.show(io :: IO, fam :: DFamily)
    ft = ftype(fam)
    if all(ft .== ft[1])
        group = "Z_$(ft[1])^$(length(fam))"
        println(io, "Family: $group -> $group")
    else
        #Z_n1 x Z_n2 x ... Z_nk
        group = ""
        for ni in ft
            group *= "Z_$(ni) × "
        end
        group = group[1 : end - 3]
        println(io, "Family: $group -> $group")
    end
end


function view(fam :: DFamily)
    ft = ftype(fam)
    k = length(fam)
    indices = CartesianIndices(ft) .- CartesianIndex(ones(Int, k)...)
    for i in indices
        x = Tuple(i)
        println(x, " -> ", fam(x...))
    end
end
