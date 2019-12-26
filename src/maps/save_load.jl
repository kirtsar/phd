using Serialization

function save_propers(ps, path_out = "../data/propers/")
	ft = ftype(ps[1])
	fname = path_out * "$ft"
	fout = open(fname, "w")
	serialize(fout, ps)
	close(fout)
end


function load_propers(ft :: Tuple, path_out = "../data/propers/")
	fname = path_out * "$ft"
	fout = open(fname, "r")
	res = deserialize(fout)
	res = Vector{typeof(res[1])}(res)
	return res
end

