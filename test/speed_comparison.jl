d = 4    # number of dimensions
K = 50   # number of points along each dimension
N = 100000  # number of points at which to interpolate

A = rand([K for i = 1:d]...)    # filtered coefficients
B = rand(N,d)                    # points at which to evaluate the splines
max(B, minimum(A)+0.01)
min(B, maximum(A)-0.01)


# splines.jl

using splines
# dummy boundaries
a = [0.0 for i=1:d]
b = [1.0 for i=1:d]
orders = Int64[K for i=1:d]
coeffs = filter_coeffs(a,b,orders,A)

# jit
values_1 = eval_UC_spline(a,b,orders,coeffs,B);
tot = sum(values_1)
tic();
values_1 = eval_UC_spline(a,b,orders,coeffs,B);
total_1 = sum(values_1)
elapsed_1 = toc();

# interpolations

using Interpolations
itp = interpolate(A, BSpline(Quadratic(Reflect())), OnGrid())



# what is the fastest way to operate on a list of points ?
LB = [copy(slice(B,i,:)) for i=1:size(B,1)]
function evaluate(itp, LB)
    s = 0.0
    for i=1:size(LB,1)
        I = LB[i]
        s += itp[I...]
    end
    return s
end

total_2 = evaluate(itp, LB)

tic();
total_2 = evaluate(itp, LB)
elapsed_2 = toc();

println("Total time: (d=",d,", K=", K, ", N=", N,")")
println("splines.jl: ",elapsed_1)
println("interpolations.jl: ",elapsed_2)
