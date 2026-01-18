module Derivatives

using ForwardDiff
# using TaylorSeries

export nth_derivative, fast_nth_derivative

# 標準版：汎用性が高く、依存が少ない
function nth_derivative(f, x, n::Int)
    n < 0 && throw(ArgumentError("n must be ≥ 0"))
    n == 0 && return f(x)

    g = f
    for _ = 1:n
        g = let current_g = g
            x -> ForwardDiff.derivative(current_g, x)
        end
    end
    return g(x)
end

# Ultrafast version
# TaylorSeries is required
function fast_nth_derivative(f, x, n::Int)
    #= 
    n < 0 && throw(ArgumentError("n must be ≥ 0"))
    n == 0 && return f(x)

    t = x + Taylor1(typeof(x), n)
    res = f(t)
    return res[n] * factorial(n)
    =#
    return nth_derivative(f, x, n)
end

end


