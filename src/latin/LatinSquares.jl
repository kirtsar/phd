module LatinSquares
	
	using AutoHashEquals
	using Combinatorics
	include("../../ProperFamilies/src/ProperFamilies.jl")
	using .ProperFamilies
	include("pairing.jl")
	include("latsquares.jl")
	include("assoc_triples.jl")
	# include("mappings.jl")

	export LatinSquare, latin, latin!
	export assoc, assoc!
	export Map

end # module
