################################
# nice printing
# helper 
# find non-zero bits in number
function nzBits(number)
    res = []
    i = 1
    while number > 0
        bit = number & 0x1
        if bit != 0
            append!(res, i)
        end
        number >>= 1
        i += 1
    end
    return res
end


# nice printing of monoms
function Base.show(io::IO, m::Monom)
    if m.val == 0
        print(io, "1")
    else
        for u in nzBits(m.val)
            print(io, "x$u")
        end
    end
end


# nice printing of zhegalkin polys
function Base.show(io::IO, f::ZhegFun)
    if isempty(f.ANF)
        print(io, "0")
    else
        for i in 1 : length(f.ANF) - 1
            mon = f.ANF[i]
            print(io, mon)
            print(io, " + ")
        end
        print(io, f.ANF[end])
    end
end


# nice printing of family of functions
function Base.show(io::IO, fam::Family)
    print(io, "Family(")
    for i in 1: length(fam.fs) - 1
        print(io, fam.fs[i])
        print(io, ", ")
    end
    print(io, fam.fs[end])
    print(io, ")")
end