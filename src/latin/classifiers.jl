using Combinatorics
using Base.Iterators   # for product
using LightGraphs
using LightGraphs.Experimental
using DataStructures  # DefaultDict


project_folder = ".."
# all possible graphs of essentials
include(project_folder * "/graphs/graph3.jl")
include(project_folder * "/graphs/graph4.jl")
# ZhegFun and Family
include("zhegalkin.jl")
include("properfam_utils.jl")

# given a family of length 3
# find the graph of its essentials
function classify3(fam :: Family)
    mat = DiGraph(graph_matrix(fam))
    for i in 1 : length(GRAPH3)
        if has_isomorph(mat, GRAPH3[i])
            return i
        end
    end
    return 0
end



function make_classes3()
    classifier = 
    DefaultDict{Int64,Vector{Family}}(() -> Vector{Family}([]))
    n, funs = generate_all_proper(3)
    for funFamily in funs
        push!(classifier[classify3(funFamily)], funFamily)
    end
    return classifier
end


# given the family, return its graph number
function classify4(fam :: Family)
    mat = DiGraph(graph_matrix(fam))
    for i in 1 : length(GRAPH4)
        if has_isomorph(mat, GRAPH4[i])
            return i
        end
    end
    return 0
end

# generate and classify all proper families of
# 4 variables
function make_classes4(funs = [])
    classifier = 
    DefaultDict{Int64,Vector{Family}}(() -> Vector{Family}([]))
    if isempty(funs)
        n, funs = experiment(4)
    end
    for funFamily in funs
        push!(classifier[classify4(funFamily)], funFamily)
    end
    return classifier
end


# given the array of families
# split by orbits
# return array of orbits
function split_orbs(fams :: Vector{Family{T}}) where T
    allFams = Set(fams)
    res = Vector([])
    while !isempty(allFams)
        fam = rand(allFams)
        orbFam = orb(fam)
        push!(res, orbFam)
        allFams = setdiff(allFams, orbFam)
    end
    return res
end


function orbs_dict(classes)
    orbsDict = Dict{Int, Array{Set{Family}}}()
    for i in keys(classes)
        orbsDict[i] = orbs(classes[i])
    end
    return orbsDict
end

#=
function print_orbs(filename, orbsDict)
    f = open(filename, "w")
    println(f, "Classification of 3-families. \n")
    for graphClass in sort(collect(keys(orbsDict)))
        println(f, "Graph number : $graphClass")
        println(f, "Orbits :")
        for orbNum in 1 : length(orbsDict[graphClass])
            print(f, "Orbit $orbNum : ")
            curr_orb = collect(orbsDict[graphClass][orbNum])
            # print current orbit 
            for fam in curr_orb
                print(f, fam)
                if fam != curr_orb[end]
                    print(f, ", ")
                end
            end
            println(f, "")
        end
        println(f, "\n")
    end
    close(f)
end
=#


function lenghts(classifier)
    m = maximum(collect(keys(classifier)))
    lens = zeros(Int, m)
    for k in keys(classifier)
        lens[k] = length(classifier[k])
    end
    return lens
end




function save_classification(clf, filename)
    out = open(project_folder * "/DATA/" * filename, "w")
    familyLen = length(clf[1][1].fs)
    nzGraphs = length(clf)

    println(out, familyLen)
    println(out, nzGraphs)
    println(out, "")

    for key in sort(collect(keys(clf)))
        # graph number
        println(out, key)
        familySet = clf[key]
        setLen = length(familySet)

        # number of families with this graph type
        println(out, setLen)
        for i in 1 : setLen
            curr_family = familySet[i]
            # current family index
            println(out, i)
            zhegArray = curr_family.fs
            for j in 1 : familyLen
                currZheg = zhegArray[j]
                currANF = currZheg.ANF
                for k in 1 : length(currANF)
                    currMonom = currANF[k]
                    val = currMonom.val
                    # print each monom in file
                    print(out, val)
                    print(out, " ")
                end
                # between different families
                println(out, "")
            end
        end
    end
    close(out)
end



function load_classification(filename)
    in = open(project_folder * "/data/" * filename, "r")
    familyLen = parse(Int, readline(in))
    nzGraphs = parse(Int, readline(in))
    # empty string
    emp = readline(in)
    classifier = 
    DefaultDict{Int64,Vector{Family{Int}}}(() -> Vector{Family{Int}}([]))

    for i in 1 : nzGraphs
        graphNum = parse(Int, readline(in))
        familysetSize = parse(Int, readline(in))
        for j in 1 : familysetSize
            familyNumber = parse(Int, readline(in))
            @assert j == familyNumber
            zhegArray = []
            for k in 1 : familyLen
                monoms = parse.(Int, split(readline(in)))
                zhegfun = ZhegFun(monoms)
                push!(zhegArray, zhegfun)
            end
            fam = Family(zhegArray...)
            push!(classifier[graphNum], fam)
        end
    end
    close(in)
    return classifier
end


function together(array_of_arrays)
    return vcat(array_of_arrays...)
end



# take one element from each orbit
function transversal(fams :: Vector{Family{T}}) where T
    famOrbs = split_orbs(fams)
    transv = Vector{Family{T}}([])
    for famOrb in famOrbs
        famElem = rand(famOrb)
        push!(transv, famElem)
    end
    return transv
end