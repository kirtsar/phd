import Base.show

using NamedArrays
using DataStructures
using StatsBase


struct Quasigroup{T}
	ls :: LatinSquare
	rowp :: Dict{Int, Perm{Int}}
	colp :: Dict{Int, Perm{Int}}
	rcnames :: Dict{T, Int}
end


# given latin square, create quasigroup 
# with this latin square as Cayley table
function Quasigroup(ls :: LatinSquare)
	n = length(ls)
	rowp, colp = mult_group(ls)
	rcnames = Dict(zip(1 : n, 0 : n - 1))
	return Quasigroup(ls, rowp, colp, rcnames)
end



function Quasigroup(ls :: LatinSquare, rcnames)
	n = length(ls)
	rowp, colp = mult_group(ls)
	#rcnames = Dict(zip(1 : n, 0 : n - 1))
	return Quasigroup(ls, rowp, colp, rcnames)
end


# given proper family and pairing
# construct quasigroup with operation
# (x , y) -> x ⊕ y ⊕ f^π(x, y)
function Quasigroup(fam :: Family, pair :: Pairing)
	ls = latin(fam, pair)
	return Quasigroup(ls)
end


# given proper family and mappings
# construct quasigroup with operation
# (x , y) -> x + y + ϕ(f(x)) + ψ(g(y))
function Quasigroup(F :: DFamily, G :: DFamily, phi :: MFamily, psi :: MFamily)
	@assert ftype(F) == ftype(G)
	@assert ftype(phi) == ftype(psi)
	@assert ftype(phi) == ftype(F)
	Gr = Grp(ftype(F))
	g_elems = elements(Gr)
	n = length(g_elems)
	rcnames = Dict(zip(g_elems, 1 : n))
	square = zeros(Int, (n, n))
	for x in g_elems
		for y in g_elems
			val = Gr(Gr(Gr(x, y), phi(F(x))), psi(G(y)))
			idx = rcnames[x]
			idy = rcnames[y]
			idv = rcnames[val]
			square[idx, idy] = idv
		end
	end
	ls = LatinSquare(square)
	res = Quasigroup(ls, rcnames)
end


function (Q :: Quasigroup)(a :: Int, b :: Int)
	return Q.ls(a, b)
end


function (Q :: Quasigroup)(a, b)
	i = rowind[a]
	j = colind[b]
	return Q(i, j)
end


function rowperm(Q :: Quasigroup, i)
	return Q.rowp[i]
end


function colperm(Q :: Quasigroup, j)
	return Q.colp[j]
end


function num_assoc(Q :: Quasigroup)
	ass_no = 0
	n = length(Q.ls)
	for i in 1 : n
		for j in 1 : n
			p1 = rowperm(Q, i)
			p2 = colperm(Q, j)
			commutator = comm(p1, p2)
			ass_no += length(fixpt(commutator))
		end
	end
	return ass_no
end


function triples(Q :: Quasigroup)
	return triples(Q.ls)
end


function Base.show(io :: IO, Q :: Quasigroup) 
	n = length(Q.ls)
	println("Quasigroup of order $n:")
	display(Q.ls.sq)
end


function full_table(Q :: Quasigroup)
	n = length(Q.ls)
	keynames = OrderedDict(Q.rcnames)
	rev = Dict([value => key for (key, value) in keynames]...)
	margin_names = [rev[key] for key in sort!(collect(keys(rev)))]
	sq = Q.ls.sq
	sqnew = [rev[key] for key in sq]
	res = NamedArray(sqnew, (margin_names, margin_names), ("A", "B"))
	return res
end



struct QuasigroupGenerator{T, S1, S2, T1, T2}
	ft :: T
	FS :: S1
	GS :: S2
	phis :: T1
	psis :: T2
end


function QuasigroupGenerator(ps)
	ft = ftype(ps[1])
	FS = ps
	GS = deepcopy(FS)
	phis = MFamilyCollector(ft)
	psis = MFamilyCollector(ft)
	return QuasigroupGenerator(ft, FS, GS, phis, psis)
end

function Base.show(io :: IO, gen :: QuasigroupGenerator)
	G = Grp(gen.ft)
	print(io, "Quasigroup Generator over $(repr(G))")
end


function Base.getindex(gen :: QuasigroupGenerator, i1, i2, j1, j2)
	F = gen.FS[i1]
	G = gen.GS[i2]
	phi = gen.phis[j1]
	psi = gen.psis[j2]
	return Quasigroup(F, G, phi, psi)
end


function ftype(gen :: QuasigroupGenerator)
	return ftype(gen.FS[1])
end


function Base.size(gen :: QuasigroupGenerator)
	return (length(gen.FS), length(gen.GS), length(gen.phis), length(gen.psis))
end


function assoc(gen :: QuasigroupGenerator)
	res = Dict{NTuple{4, Int}, Int}()
	inds = CartesianIndices(size(gen))
	for ind in tqdm(inds)
		Q = gen[Tuple(ind)...]
		n_as = num_assoc(Q)
		res[Tuple(ind)] = n_as
	end
	return res
end


struct AssocAnalysis{T, S}
	ks :: T
	vals :: S
end

function getindex(a :: AssocAnalysis, i :: Int)
    return (a.ks[i], a.vals[i])
end


function getindex(a :: AssocAnalysis, I) 
    return [a[i] for i in I]
end


function minimal(a :: AssocAnalysis)
    i = argmin(a.vals)
    return a[i]
end


function minimal(a :: AssocAnalysis, i :: Int)
    inds = sortperm(a.vals)[1 : i]
    return a[inds]
end


function AssocAnalysis(d :: Dict)
	ks = collect(keys(d))
	vals = [d[key] for key in ks]
	return AssocAnalysis(ks, vals)
end


function minimal_num(as :: AssocAnalysis)
	m = minimum(as.vals)
	return count(as.vals .== m)
end



function summary(as :: AssocAnalysis)
	return countmap(as.vals)
end


# given numbers of discrete families
# find number of assoc triples
function summary(as :: AssocAnalysis, i1, i2)
	res = Int[]
	for (i, key) in enumerate(as.ks)
		if (key[1] == i1) && (key[2] == i2)
			push!(res, as.vals[i])
		end
	end
	return countmap(res)
end


function meanvalue(t :: Dict)
	res = 0.0
	total = 0
	for k in keys(t)
		res += k * t[k]
		total += t[k]
	end
	return res/total
end


function mean_table(as :: AssocAnalysis, proper_num)
	res = zeros((proper_num, proper_num))
	for i in 1 : proper_num
		for j in 1 : proper_num
			v = summary(as, i, j)
			res[i, j] = meanvalue(v)
		end
	end
	return res 
end


