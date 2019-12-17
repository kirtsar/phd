struct Grp{T}
	ftype :: T
end


function (G :: Grp)(xs, ys)
	n = length(xs)
	res = zeros(Int, n)
	for i in 1 : n
		res[i] = (xs[i] + ys[i]) % G.ftype[i]
	end
	return Tuple(res)
end


function elements(G :: Grp{T}) where T
	ft = G.ftype
	return Tuple.(vec(CartesianIndices(ft) .- CartesianIndex(ones(Int, length(ft))...)))
end


function show(io :: IO, G :: Grp)
	ft = G.ftype
	maptype = ""
	for ni in ft
		maptype *= "Z_$(ni) × "
	end
	print(io, "$(maptype[1 : end - 3])")
end