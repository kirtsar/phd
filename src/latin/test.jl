include("classifiers.jl")
clf3 = load_classification("/DATA/classification3.clf")
fam = clf3[16][1]
pi = Pairing(fst)
as_perm(fam, pi, 1)


function test_ortho()
    #=
    Family(
    ZhegFun([[3,4], [2,3,4]]),
    ZhegFun([[1,3,4], [1, 4]]),
    ZhegFun([[1,2], [1,2,4]]),
    ZhegFun([[2,3], [1,2,3]]),
    )
    =#

    fam2 = Family(
    ZhegFun([[3,4,5], [2,3,4,5], Int[]]),
    ZhegFun([[4,5,1], [3,4,5,1]]),
    ZhegFun([[5,1,2], [4,5,1,2]]),
    ZhegFun([[1,2,3], [5,1,2,3]]),
    ZhegFun([[2,3,4], [1,2,3,4]])
    )
    @assert(is_ortho(fam2) == false)
    @assert(has_ortho(fam2) == true)
    return true
end

