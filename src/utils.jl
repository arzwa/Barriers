pairfst_unidirectional_island = eval(
    :(function (m, N_A, N_B)
      #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =# @inbounds begin
              #= /home/arzwa/.julia/packages/Symbolics/ociKG/src/build_function.jl:124 =#
              begin
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:389 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:390 =#
                  #= /home/arzwa/.julia/packages/SymbolicUtils/6fncq/src/code.jl:391 =#
                  (/)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(1//2, N_A), m)), (*)((*)(1//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-1//1, N_A), (^)(m, 2))), (*)((*)(-1//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)(1//2, N_A), (^)(m, 3))), (*)((*)(1//2, N_B), (^)(m, 3))), (+)((+)((+)((+)((+)((+)((+)((+)((+)((+)((+)(1//1, (*)(-3//1, m)), (*)((*)(3//2, N_A), m)), (*)((*)(7//2, N_B), m)), (*)(3//1, (^)(m, 2))), (*)((*)(-3//1, N_A), (^)(m, 2))), (*)((*)(-5//1, N_B), (^)(m, 2))), (*)(-1//1, (^)(m, 3))), (*)((*)((*)(4//1, N_A), N_B), (^)(m, 2))), (*)((*)(3//2, N_A), (^)(m, 3))), (*)((*)(3//2, N_B), (^)(m, 3))), (*)((*)((*)(-2//1, N_A), N_B), (^)(m, 3))))
              end
          end
  end)
)

function winstat(stat, winsize, xs, ys)
    T = typeof(stat(ys))
    zs = []
    ywin = T[]
    win1 = winsize
    i = 1
    while i <= length(xs)
        while i <= length(xs) && xs[i] < win1
            push!(ywin, ys[i])
            i += 1
        end
        win = ((win1-winsize), win1)
        push!(zs, (win, stat(ywin)))
        win1 += winsize
        ywin = T[]
    end
    return zs
end

"""
    windowed(data::AbstractVector{Union{AbstractVector,Tuple}}, windows)

The first entry of each data point in `data` shoudl be the location, and
`windows` corresponds to window endpoints. This returns the data in `data` in
bins.
"""
function windowed(data::Vector{T}, windows) where T
    i = 1
    windata = [T[] for _=1:length(windows)]
    for (x, y) in data
        x > windows[i] && (i += 1)
        push!(windata[i], (x,y))
    end
    return windata
end

