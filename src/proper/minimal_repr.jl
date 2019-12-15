function nzLen(zheg :: ZhegFun{T}) where T
    reslen = 0
    for mon in zheg.ANF
        reslen += length(nzBits(mon.val))
    end
    return reslen
end


function nzLen(fam :: Family{T}) where T
    n = length(fam)
    reslen = zeros(Int, n)
    for i in 1 : n
        reslen[i] = nzLen(fam[i])
    end 
    return reslen
end


function nz_summary(zheg :: ZhegFun{T}) where T
    reslen = 0
    for mon in zheg.ANF
        reslen += mon.val
    end
    return reslen
end


function nz_summary(fam :: Family{T}) where T
    n = length(fam)
    reslen = zeros(Int, n)
    for i in 1 : n
        reslen[i] = nz_summary(fam[i])
    end 
    return reslen
end


function minimal_repr(fam)
    n = length(fam)
    equiv_fams = collect(orb(fam))
    minfam = 0
    minlens = 2^n * ones(Int, n)
    for fam in equiv_fams
        currlens = nzLen(fam)
        if currlens < minlens
            minfam = fam
            minlens = currlens
        elseif currlens == minlens
            currSummary = nz_summary(fam)
            minSummary = nz_summary(minfam)
            if currSummary < minSummary
                minfam = fam
            end
        end
    end
    return minfam
end