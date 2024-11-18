const states = ["F", "FD", "HA", "HB", "HAB"]

function classify(x)
    (x[1] == x[2] == x[3] == x[4]) && return "F"
	(x[1] == x[2] && x[3] == x[4]) && return "FD"
	(x[1] != x[2] && x[3] == x[4]) && return "HA"
	(x[1] == x[2] && x[3] != x[4]) && return "HB"
	(x[1] != x[2] && x[3] != x[4]) && return "HAB"
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
