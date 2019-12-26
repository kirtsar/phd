import Base.*
import Base.inv
import Base.show
import Base.getindex
import Base.setindex!

using OffsetArrays
using Base.Iterators
using IterTools


abstract type AbstractMap end

# mapping Zk -> Zk
# represent as array
# f(k) = img[k]
struct Map{T} <: AbstractMap
	img :: OffsetVector{T}
	ft :: Int
end


# bijective map {1, .. k} -> {1, ... k}
struct Perm{T} <: AbstractMap
	img :: Vector{T}
	ft :: Int
end

# mapping of the form {0, .. k - 1} -> {0, .., k - 1}
# NO  CHECKS
function Map(x :: AbstractArray)
	v = OffsetVector(x, 0 : length(x) - 1)
	return Map(v, length(x))
end


function Map(x :: Tuple)
	xarr = collect(x)
	return Map(xarr)
end


function Perm(x :: AbstractArray)
	ft = length(x)
	return Perm(x, ft)
end


function Map(f :: Family)
	n = length(f)
	x = zeros(Int, 2^n)
	for i in 1 : 2^n
		x[i] = f(i - 1)
	end
	return Map(x)
end


function (m :: AbstractMap)(x)
	return m.img[x]
end


function ftype(m :: AbstractMap)
	return m.ft
end


function setindex!(m :: AbstractMap, val, i :: Int)
	m.img[i] = val
end


function Base.show(io :: IO, m :: Map)
	ft = ftype(m)
	print(io, "Map Z_$ft -> Z_$ft")
end


function Base.show(io :: IO, m :: Perm)
	ft = ftype(m)
	print(io, "Permutation on {1, ..., $ft}")
end


function view(m :: AbstractMap)
	for i in eachindex(m.img)
		println(i, " -> ", m(i))
	end
end

#=
function view(m :: Perm)
	ft = ftype(m)
	for i in 1 : ft
		println(i, " -> ", m(i))
	end
end
=#

function *(m1 :: Perm{T}, m2 :: Perm{T}) where T
	n = ftype(m1)
	x = zeros(T, n)
	for i in 1 : n
		# i -> m2(i) -> m1(m2(i))
		x[i] = m1.img[m2.img[i]]
	end
	return Perm(x)
end


function inv(m :: Perm)
	n = ftype(m)
	x = zeros(Int, n)
	for i in 1 : n
		x[m(i)] = i
	end
	return Perm(x)
end


function fixpt(m :: Perm{T}) where T
	fix = T[]
	for i in eachindex(m.img)
		if m(i) == i
			push!(fix, i)
		end
	end
	return fix
end


function comm(m1 :: Perm, m2 :: Perm)
	return m1 * m2 * inv(m1) * inv(m2)
end


function mult_group(ls :: LatinSquare)
	rowperm = Dict{Int, Perm{Int}}()
	colperm = Dict{Int, Perm{Int}}()
	n = length(ls)
	for i in 1 : n
		rowperm[i] = Perm(ls.sq[i, :])
		colperm[i] = Perm(ls.sq[:, i])
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


struct MapCollector{T}
	ms :: Vector{Map{T}}
end


function MapCollector(n :: Int)
	inds = CartesianIndices(Tuple(n * ones(Int, n))) .- CartesianIndex(Tuple(ones(Int, n)))
	vecs = vec(Tuple.(inds))
	ms = Map.(vecs)
	return MapCollector(ms)
end


function Base.show(io :: IO, mcol :: MapCollector)
	print(io, "Collection of Maps: $(repr(mcol.ms[1]))")
end


struct MFamilyCollector{T}
	coll :: T
end


function MFamilyCollector(ft :: Tuple)
	n = length(ft)
	coll = Any[undef for i in 1 : n]
	for i in 1 : n 
		coll[i] = MapCollector(ft[i]).ms
	end
	mfs = Iterators.product(coll...)
	return MFamilyCollector(mfs)
end


function Base.getindex(mcol :: MFamilyCollector, n)
	return MFamily(collect(nth(mcol.coll, n)))
end


function Base.length(mcol :: MFamilyCollector)
	return length(mcol.coll)
end



using ProgressBars
# ps = vector of proper families
function generate_all_quasigroups(ps)
	res = []
	ft = ftype(ps[1])
	mcols1 = MFamilyCollector(ft)
	mcols2 = deepcopy(mcols1)
	n = length(ps)
	for i in tqdm(1 : n)
		for j in i : n
			F = ps[i]
			G = ps[j]
			for mcol1 in mcols1.coll
				for mcol2 in mcols2.coll
					phi = MFamily(collect(mcol1))
					psi = MFamily(collect(mcol2))
					Q = Quasigroup(F, G, phi, psi)
					push!(res, Q)
				end
			end
		end
	end
	return res
end



#=

=#





