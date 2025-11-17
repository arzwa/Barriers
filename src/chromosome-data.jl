"""
    ChromosomeData

All the data we need to conduct an analysis.

- snps : A DataFrame containing biallelic SNPs in rows
- acc  : a vector of accessible intervals 
- len  : Total chromosome length
- popA : indices/colnames of population A (mainland)
- popB : indices/colnames of population B (island)
- xcol : index/colname of the SNP positions

"""
struct ChromosomeData{T<:AbstractVector,V,W}
    snps :: DataFrame
    acc  :: Vector{Tuple{W,W}}
    len  :: W
    popA :: T
    popB :: T
    xcol :: V
end

"""
    isoverlapping(int1, int2)

Find out whether two **closed** intervals are overlapping, i.e. have at least
one element in common.
"""
function isoverlapping(int1, int2)
    x0, x1 = extrema(int1)
    y0, y1 = extrema(int2)
    x0 <= y1 && y0 <= x1
end

"""
    overlaplen(int1, int2)

Calculate overlap length of two intervals.
"""
function overlaplen(int1, int2)
    x0, x1 = extrema(int1)
    y0, y1 = extrema(int2)
    !(x0 <= y1 && y0 <= x1) && return 0
    (x0 >= y0 && x1 <= y1)  && return length(int1)
    (x0 <= y0 && x1 >= y1)  && return length(int2)
    (x0 >= y0 && x1 >= y1)  && return y1 - x0 + 1
    (x0 <= y0 && x1 <= y1)  && return x1 - y0 + 1
end

get_accsites_windows(chrom::ChromosomeData, windows) = 
    get_accsites_windows(chrom.acc, windows)
function get_accsites_windows(acc, windows)
    i = 1  # index going over accessible stretches
    map(windows) do (w0, w1)
        n = 0  # counter for number of accessible sites in the window
        # note that windows are closed intervals (w0 and w1 are included in the window)
        while i <= length(acc) && (acc[i][end] < w0 || isoverlapping(acc[i], w0:w1))
            n += overlaplen(acc[i], w0:w1)
            acc[i][end] > w1 && break  # if the acc stretch extends beyond the windo, don't advance!
            i += 1
        end
        n
    end
end

function _numcombos(chrom::ChromosomeData)
    @unpack popA, popB = chrom
    a = length(popA); b = length(popB)
    a*(a-1)*b*(b-1)÷4
end

"""
    snpstats

Calculate per-site statistics for the two populations (allele frequencies, pi
within, pi between, fst).
"""
function snpstats(data::ChromosomeData)
    @unpack snps, popA, popB, xcol = data
    sdf = select(snps, xcol, AsTable(:) => ByRow(x->_getsitestats(x, popA, popB)) => AsTable)
end

function _getsitestats(x, A, B)
    as = [x[a] for a in A if !ismissing(x[a])]
    bs = [x[b] for b in B if !ismissing(x[b])]
    j = sum(as); m=length(as)
    k = sum(bs); n=length(bs)
    pa  = j/m  # all. freq. A
    pb  = k/n  # all. freq. B
    a01 = j*(m-j)
    a00 = j*(j-1)÷2 
    a11 = (m-j)*(m-j-1)÷2
    b01 = k*(n-k)
    b00 = k*(k-1)÷2
    b11 = (n-k)*(n-k-1)÷2
    T   = m*(m-1)*n*(n-1)÷4
    HAB = a01*b01
    HA  = a01*(b00+b11)
    HB  = b01*(a00+a11)
    FD  = a00*b11 + a11*b00
    F   = T - HAB - HA - HB - FD
    piA = a01/(a00 + a01 + a11)
    piB = b01/(b00 + b01 + b11)
    pib = (j/m)*(n-k)/n + (m-j)/m * k/n
    piw = (length(A)*piA + length(B)*piB)/(length(A) + length(B))
    fst = (pib - piw)/(pib+piw)
    (piA=piA, piB=piB, piw=piw, pib=pib, fst=fst,
        F=F, FD=FD, HA=HA, HB=HB, HAB=HAB)
end

function windowdata(chrom::ChromosomeData, windows, df=snpstats(chrom))
    @unpack len, xcol = chrom
    ns = get_accsites_windows(chrom, windows)
    nc = _numcombos(chrom)
    xdf = Barriers.window_df(df, windows, fun=x->sum(filter(!isnan, skipmissing(x))), xcol=xcol)
    xdf = hcat(xdf, DataFrame(:ns=>ns))
    sstats = [:piA, :piB, :piw, :pib, :fst]
    cstats = [:F, :FD, :HA, :HB, :HAB]
    _fun1(x) = (; [k=>x[k]/x[:ns] for k in sstats]...)
    _fun2(x) = (x[:F] + nc*x[:ns] - sum([x[k] for k in cstats]))
    select!(xdf, 
        chrom.xcol, 
        #:w0, 
        :w1, :wmid, :winlen, :ns,
        AsTable([:ns; sstats]) => ByRow(_fun1) => AsTable,
        AsTable([:ns; cstats]) => ByRow(_fun2) => :F,
        :FD, :HA, :HB, :HAB
    )
    transform!(xdf, 
        AsTable([:piw,:pib]) => ByRow(x->(x[:pib]-x[:piw])/(x[:pib]+x[:piw])) => :fstw
    )
    return xdf
end

getwindows(chrom::ChromosomeData, ws) = getwindows(chrom.len, ws)
function getwindows(len::Int, ws)
    xs = [collect(1:ws:len) ; len+1]
    [(xs[i], xs[i+1]-1) for i=1:length(xs)-1]
end

function getwindows(len, ws)
    n = len ÷ ws + 1
    [((i-1)*ws, min(i*ws, len)) for i=1:n]
end

"""
    window_df

Transform a DataFrame with site statistics to one with window statistics, using
`fun` to aggregate the site values within windows.
"""
function window_df(df, windows; fun=x->mean(filter(!isnan, skipmissing(x))), xcol=:POS)
    xs = map(eachcol(df)) do col
        _winstat(fun, windows, df[:,xcol], col)
    end
    wdf = DataFrame(hcat(xs...), names(df))
    x0 = first.(windows)
    x1 = last.(windows)
    wdf = hcat(DataFrame(:w0=>x0, :w1=>x1, :wmid=>mean.(windows), :winlen=>x1 .- x0 .+ 1), wdf) 
    return wdf
end

function _winstat(stat, windows, xs, ys::AbstractVector{T}) where T
    V = typeof(stat(ys))
    zs = Vector{V}(undef, length(windows))
    ywin = T[]
    i = 1  # site index
    j = 1  # window index
    while j <= length(windows)
        win0, win1 = windows[j]
        while i <= length(xs) && xs[i] <= win1
            push!(ywin, ys[i])
            i += 1
        end
        zs[j] = stat(ywin)
        j += 1
        ywin = T[]
    end
    return zs
end




