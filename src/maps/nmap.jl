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


function (nm :: NMap)(inds...)
	return nm[inds...]
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

#=
function Base.*(m :: Map, nm :: NMap)
	new_img = deepcopy(nm.img)
=#


