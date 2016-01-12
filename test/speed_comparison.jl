d = 3    # number of dimensions
K = 50   # number of points along each dimension
N = 100000  # number of points at which to interpolate

A = rand([K for i = 1:d]...)    # filtered coefficients
B = rand(N,d)                    # points at which to evaluate the splines
max(B, minimum(A)+0.01)
min(B, maximum(A)-0.01)

println("Comparison (d=",d,", K=", K, ", N=", N,")")


# interpolations.jl
using Interpolations
itp = interpolate(A, BSpline(Quadratic(Reflect())), OnGrid())

# what is the fastest way to operate on a list of points ?
LB = Vector{Float64}[copy(slice(B,i,:)) for i=1:size(B,1)]




function evaluate(itp, LB)
    s = 0.0
    for i=1:size(LB,1)
        I = LB[i]
        s += itp[I[1],I[2],I[3]]
    end
    return s
end

total_2 = evaluate(itp, LB)

println("interpolations.jl (quadratic)")
@time total_2 = evaluate(itp, LB)
exit()


# splines.jl
using splines

# dummy boundaries
a = [0.0 for i=1:d]
b = [1.0 for i=1:d]
orders = Int64[K for i=1:d]
coeffs = filter_coeffs(a,b,orders,A)

# jit
println("splines.jl: ")
values_1 = eval_UC_spline(a,b,orders,coeffs,B);
tot = sum(values_1)
@time  values_1 = eval_UC_spline(a,b,orders,coeffs,B); tot = sum(values_1)
# not sure whether the summation is taken into accoun
println("splines.jl (with derivatives): ")
values_1, dvalues_1 = eval_UC_spline_G(a,b,orders,coeffs,B); tot = sum(values_1)

@time  values_1, dvalues_1 = eval_UC_spline_G(a,b,orders,coeffs,B);total_1 = sum(values_1)
