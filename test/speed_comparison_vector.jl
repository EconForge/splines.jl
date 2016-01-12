d = 3    # number of dimensions
K = 50   # number of points along each dimension
N = 10^6  # number of points at which to interpolate
vecsize = 10

orders = Int64[K for i=1:d]
A = rand([K for i = 1:d]...)    # filtered coefficients
B = rand(N,d)                    # points at which to evaluate the splines


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

println("interpolations.jl (cubic)")
@time total_2 = evaluate(itp, LB)
# exit()


using FixedSizeArrays

AA = rand(vecsize,orders...)
A0 = reinterpret(Vec{vecsize,Float64}, AA, (K,K,K))

bitp = interpolate(A0, BSpline(Cubic(Flat())), OnGrid())
function foo(itp, X)
    s = zero(eltype(itp))
    for x in X
        # s += itp[x[1]]
        s += itp[x[1],x[2],x[3]]
    end
    s
end

foo(bitp, LB)
println("interpolations.jl: (vector x",vecsize,")")

@time foo(bitp, LB)
# exit(0)
#
# orders = Int64[K for i in 1:d]
# a = rand(vecsize,orders...)
# b = reinterpret(Vec{vecsize,Float64}, a, (K,K,K))
# bitp = interpolate(b, BSpline(Cubic(Flat())), OnGrid())
# function foo(itp, X)
#     s = zero(eltype(itp))
#     for x in X
#         # s += itp[x[1]]
#         s += itp[x[1],x[2],x[3]]
#     end
#     s
# end
#
# # X = 2+150*rand(N)
# foo(bitp, XX) # warmup
# println("interpolations.jl (vector x",vecsize,")")
# @time foo(bitp, XX)


##################skyp
# with splines.jl #
###################


# splines.jl
using splines

# dummy boundaries
a = [0.0 for i=1:d]
b = [1.0 for i=1:d]

vals = rand(orders...)
coeffs = filter_coeffs(a,b,orders,vals)

# mcoeffs represent vecsize identical splines
mcoeffs = zeros(vecsize, size(coeffs)...)
for s in 1:vecsize
    mcoeffs[s,:] =  coeffs
end
# jit
# X = 150*rand(N,d)


println("splines.jl: (scalar)")
values_1 = eval_UC_spline(a,b,orders,coeffs,B);
tot = sum(values_1)
@time  values_1 = eval_UC_spline(a,b,orders,coeffs,B); tot = sum(values_1)

println("splines.jl: (vector x",vecsize,")")
values_1 = eval_UC_multi_spline(a,b,orders,mcoeffs,B);
tot = sum(values_1)
@time  values_1 = eval_UC_multi_spline(a,b,orders,mcoeffs,B); tot = sum(values_1)
@time  eval_UC_multi_spline!(a,b,orders,mcoeffs,B, values_1);
