import Base.show

struct Quasigroup
	ls :: LatinSquare
	rowp :: Dict{Int, Map{Int}}
	colp :: Dict{Int, Map{Int}}
end


function Quasigroup(ls :: LatinSquare)
	rowp, colp = mult_group(ls)
	return Quasigroup(ls, rowp, colp)
end


function Quasigroup(fam :: Family, pair :: Pairing)
	ls = latin(fam, pair)
	return Quasigroup(ls)
end


function (Q :: Quasigroup)(a, b)
	return Q.ls(a, b)
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


function Base.show(io :: IO, Q :: Quasigroup) 
	n = length(Q.ls)
	println("Quasigroup of order $n:")
	display(Q.ls.sq)
end