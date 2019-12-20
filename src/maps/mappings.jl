import Base.*
import Base.inv
import Base.show
import Base.setindex!

using OffsetArrays


# mapping Zk -> Zk
# represent as array
# f(k) = img[k]
struct Map{T}
	img :: OffsetVector{T}
	k :: Int
end


# mapping of the form {0, .. k - 1} -> {0, .., k - 1}
# NO  CHECKS
function Map(x :: AbstractArray)
	v = OffsetVector(x, 0 : length(x) - 1)
	return Map(v, length(x))
end


function Map(f :: Family)
	n = length(f)
	x = zeros(Int, 2^n)
	for i in 1 : 2^n
		x[i] = f(i - 1)
	end
	return Map(x)
end


function (m :: Map)(x)
	return m.img[x]
end


function ftype(m :: Map)
	return m.k
end


function setindex!(m :: Map, val, i :: Int)
	m.img[i] = val
end


function Base.show(io :: IO, m :: Map)
	ft = ftype(m)
	print(io, "Map Z_$ft -> Z_$ft")
end


function view(m :: Map)
	ft = ftype(m)
	for i in 0 : ft - 1
		println(i, " -> ", m(i))
	end
end


function *(m1 :: Map{T}, m2 :: Map{T}) where T
	n = ftype(m1)
	x = zeros(T, n)
	for i in 0 : n-1
		# i -> m2(i) -> m1(m2(i))
		x[i+1] = m1.img[m2.img[i]]
	end
	return Map(x)
end


function inv(m :: Map)
	n = ftype(m)
	x = zeros(Int, n)
	for i in 0 : n - 1
		x[m(i) + 1] = i
	end
	return Map(x)
end


function fixpt(m :: Map{T}) where T
	n = len(m)
	fix = T[]
	for i in 0 : n - 1
		if m(i) == i
			push!(fix, i)
		end
	end
	return fix
end


function comm(m1 :: Map, m2 :: Map)
	return m1 * m2 * inv(m1) * inv(m2)
end


function mult_group(ls :: LatinSquare)
	rowperm = Dict{Int, Map{Int}}()
	colperm = Dict{Int, Map{Int}}()
	n = length(ls)
	for i in 1 : n
		rowperm[i - 1] = Map(ls.sq[i, :])
		colperm[i - 1] = Map(ls.sq[:, i])
	end
	return rowperm, colperm
end



# collection of mappings
# f1 : Zk1 -> Zk1
# ...
# ft : Zkt -> Zkt
struct MFamily{T, S}
	fs :: Vector{Map{T}}
	ftype :: S
end


function MFamily(ms :: Vector{Map{T}}) where T
	ft = Tuple(ftype.(ms))
	return MFamily(ms, ft)
end


function getindex(mfam :: MFamily, i :: Int)
	return mfam.fs[i]
end


function (mfam :: MFamily)(x :: Tuple)
	t = length(x)
	res = zeros(Int, t)
	# f1(x1), ..., ft(xt)
	for i in 1 : t
		res[i] = mfam[i](x[i])
	end
	return Tuple(res)
end


function ftype(mfam :: MFamily)
	return mfam.ftype
end


function Base.show(io :: IO, m :: MFamily)
	ft = ftype(m)
	G = Grp(ft)
	maptype = repr(G)
	print(io, "Map $maptype -> $maptype")
end


function view(m :: MFamily)
	ft = ftype(m)
	G = Grp(ft)
	xs = elements(G)
	for x in xs
		println(x, " -> ", m(x))
	end
end

