import Base.*
import Base.inv
import Base.show

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


function Base.show(io :: IO, nm :: NMap)
	n = size(nm.img)
	if all(n .== n[1])
		println("Map of the type Z_$(n[1])^$(length(n))")
	else
		#Z_n1 x Z_n2 x ... Z_nk
		maptype = ""
		for ni in n
			maptype *= "Z_$(ni) × "
		end
		println("Map of the type $(maptype[1 : end - 3])")
	end
end


function (nm :: NMap)(inds :: NTuple)
	return nm.img[inds...]
end


function view(nm :: NMap)
	n = size(nm.img)
	indices = CartesianIndices(n) .- CartesianIndex(1, 1, 1)
	for i in indices
		xi = Tuple(i)
		println(xi, " : ", nm(xi))
	end
end

#=
function Base.*(m :: Map, nm :: NMap)
	new_img = deepcopy(nm.img)
=#


