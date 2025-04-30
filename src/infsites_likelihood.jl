
# Migration is A->B forward in time
function CLModel(m, μ, a, b, n)
    p = probs(m, μ, a, b)
    return Multinomial(n, p)
end

const states = ["F", "FD", "HA", "HB", "HAB"]

#function logpdfcl(m, μ, a, b, x)
#    p = probs(m, μ, a, b)
#    if !isprobvec(p) 
#        @warn "Not a probability vector!"
#        return -Inf
#    else
#        return logpdf(Multinomial(sum(x), p), x)
#    end
#end

logpdfsite(m, u, a, b, site) = log(probs(m, u, a, b)[site])

# γ is the number of pseudo-data points per window
function _logpdfcl(m, μ, a, b, x, γ=1.0)
    p = siteprobs(m, μ, a, b)
    y = x .* γ  # pseudo counts
    n = sum(y)
    loggamma(sum(y) + 1) + sum((y .* log.(p))) - sum(loggamma.(y .+ 1))
end  # Multinomial

function logpdfcl(m, μ, a, b, x, γ=1.0)
    p = siteprobs(m, μ, a, b)
    y = x .* γ  # pseudo counts
    sum(log.(p) .* y)
end  # Categorical

function classify(x)
    (x[1] == x[2] == x[3] == x[4]) && return "F"
	(x[1] == x[2] && x[3] == x[4]) && return "FD"
	(x[1] != x[2] && x[3] == x[4]) && return "HA"
	(x[1] == x[2] && x[3] != x[4]) && return "HB"
	(x[1] != x[2] && x[3] != x[4]) && return "HAB"
end

classify2(x) = findfirst(y->y==classify(x), states)

function getcounts(L, x)
    cts = countmap(x)
    y = Dict(k => 0 for k in states)
    for (k,v) in cts
        x = classify(k)
        y[x] += v
    end
    y["F"] += L - sum(values(y))
    [y[k] for k in states]
end

function get_counts_windows(data::Vector{<:Tuple}, winsize)
    wins = Vector{Int64}[]
    curr = winsize
    temp = String[]
    i = 1
    while i <= length(data)
        if data[i][1] > curr
            cm = countmap(temp)
            for k in states
                !(haskey(cm, k)) && (cm[k] = 0)
            end
            cm["F"] += winsize - sum(values(cm))
            push!(wins, [cm[k] for k in states])
            temp = String[]
            curr += winsize
        end
        push!(temp, data[i][2])
        i += 1
    end
    cm = countmap(temp)
    for k in states
        !(haskey(cm, k)) && (cm[k] = 0)
    end
    cm["F"] += winsize - sum(values(cm))
    push!(wins, [cm[k] for k in states])
    return wins
end

# this one assumes that fixed sites are included (e.g. from Fwd sim)
function get_counts_windows2(data::Vector{<:Tuple}, winsize)
    wins = Vector{Int64}[]
    curr = winsize
    temp = String[]
    i = 1
    while i <= length(data)
        if data[i][1] > curr
            cm = countmap(temp)
            for k in states
                !(haskey(cm, k)) && (cm[k] = 0)
            end
            push!(wins, [cm[k] for k in states])
            temp = String[]
            curr += winsize
        end
        push!(temp, data[i][2])
        i += 1
    end
    cm = countmap(temp)
    for k in states
        !(haskey(cm, k)) && (cm[k] = 0)
    end
    push!(wins, [cm[k] for k in states])
    return wins
end


sitepr(m, μ, λa, λb, x=[["a","a"],["b","b"]]) = 
    sitepr(promote(m, μ, λa, λb)..., x=[["a","a"],["b","b"]])

function sitepr(m::T, μ::T, λa::T, λb::T, x=[["a","a"],["b","b"]]) where T
    res = Dict{Vector{Vector{String}},T}()
	function recursion(x, p)
		xa, xb = x
		na, nb = length.(x)
		ka = na*(na-1)//2
		kb = nb*(nb-1)//2
		rates = [na*μ, nb*μ, nb*m, λa*ka, λb*kb]
		probs = rates ./ sum(rates)
		# mutation in pop A
		for i=1:na
			_xa = copy(xa) 
			_xa[i] = uppercase(_xa[i])
			_x = [sort(_xa), sort(xb)]
			if haskey(res, _x) 
				res[_x] += p*probs[1]/na
			else
				res[_x] = p*probs[1]/na
			end
		end
		# mutation in pop B
		for i=1:nb
			_xb = copy(xb) 
			_xb[i] = uppercase(_xb[i])
			_x = [sort(xa), sort(_xb)]
			if haskey(res, _x) 
				res[_x] += p*probs[2]/nb
			else
				res[_x] = p*probs[2]/nb
			end
		end
		# migration B -> A (backwards in time)
		for i=1:nb
			_xa = copy(xa)
			_xb = copy(xb)
			mig = popat!(_xb, i)
			push!(_xa, mig)
			recursion([_xa, _xb], p*probs[3]/nb)
		end
		# coalescence in A
		for i=1:na-1
			for j=i+1:na
				_xa = copy(xa)
				xj = popat!(_xa, j)
				_xa[i] = join(sort([_xa[i], xj]), "")
				recursion([_xa, xb], p*probs[4]/ka)
			end
		end
		# coalescence in B
		for i=1:nb-1
			for j=i+1:nb
				_xb = copy(xb)
				xj = popat!(_xb, j)
				_xb[i] = join(sort([_xb[i], xj]), "")
				recursion([xa, _xb], p*probs[5]/kb)
			end
		end
	end
	recursion(x, 1)
    d = Dict{String,T}(x=>zero(T) for x in states)
	for (k,v) in res
		x = _classify(k)
        d[x] = d[x] + v
	end
	return d 
end

function _classify(x)
	xx = join(sort(vcat(x...)), "")
	A = occursin("A",xx)
	a = occursin("a",xx)
	B = occursin("B",xx)
	b = occursin("b",xx)
	((A && a) && (B && b)) && return "HAB"
	(A && a) && return "HA"
	(B && b) && return "HB"
	F = (!(A && a) && !(B && b))
	(F && ((A && b) || (a && B))) && return "FD"
	return "F"
end
