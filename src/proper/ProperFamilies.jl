module ProperFamilies

	include("zhegalkin.jl")
	include("family.jl")
	include("output_show.jl")
	include("proper.jl")
	include("group_action.jl")
	include("orthogonality.jl")
	include("minimal_repr.jl")

	export Monom, ZhegFun, Family
	export monom
	export permute, orb, stab, additive_orb
	export is_proper
	export is_ortho, has_ortho, ortho_example
	export minimal_repr

end # module
