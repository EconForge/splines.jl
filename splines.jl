module splines

    export filter_coeffs, interpolant_cspline
    export eval_UC_spline, eval_UC_spline_G
    include("csplines.jl")
    include("splines_filter.jl")

function interpolant_cspline(smin, smax, orders, V)

    coefs = filter_coeffs(smin, smax, orders, V)

    function fun(s::Array{Float64,2})
        return eval_UC_spline(smin, smax, orders, coefs, s)
    end

    function fun(p::Float64...)
        return fun([p...]')
    end
    
    return fun 

end



end
