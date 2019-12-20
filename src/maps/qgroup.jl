import Base.show

using NamedArrays
using DataStructures


struct Quasigroup{T}
	ls :: LatinSquare
	rowp :: Dict{Int, Map{Int}}
	colp :: Dict{Int, Map{Int}}
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
	for i in 0 : (n - 1)
		for j in 0 : (n - 1)
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

