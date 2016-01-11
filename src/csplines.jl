include("macros.jl")

A44d = [
   [-1.0/6.0  3.0/6.0 -3.0/6.0 1.0/6.0];
   [ 3.0/6.0 -6.0/6.0  0.0/6.0 4.0/6.0];
   [-3.0/6.0  3.0/6.0  3.0/6.0 1.0/6.0];
   [ 1.0/6.0  0.0/6.0  0.0/6.0 0.0/6.0]
]

dA44d = [
   [ 0.0 -0.5  1.0 -0.5];
   [ 0.0  1.5 -2.0  0.0];
   [ 0.0 -1.5  1.0  0.5];
   [ 0.0  0.5  0.0  0.0]
]

d2A44d = [
   [ 0.0 0.0 -1.0  1.0];
   [ 0.0 0.0  3.0 -2.0];
   [ 0.0 0.0 -3.0  1.0];
   [ 0.0 0.0  1.0  0.0]
]

function eval_UC_spline(smin, smax, orders, C, S)

    d = size(S,2)
    N = size(S,1)

    vals = zeros(N)

    if d == 1
        eval_UC_spline_1d(smin, smax, orders, C, S, vals, A44d, dA44d)
    elseif d == 2
        eval_UC_spline_2d(smin, smax, orders, C, S, vals, A44d, dA44d)
    elseif d == 3
        eval_UC_spline_3d(smin, smax, orders, C, S, vals, A44d, dA44d)
    elseif d == 4
        eval_UC_spline_4d(smin, smax, orders, C, S, vals, A44d, dA44d)
    end

    return vals

end

function eval_UC_spline_G(a, b, orders, C, S)

    d = size(S,2)
    N = size(S,1)

    vals = zeros(N)
    grad = zeros(N,d)

    if d == 1
        eval_UC_spline_1d(a, b, orders, C, S, vals, grad, A44d, dA44d)
    elseif d == 2
        eval_UC_spline_2d(a, b, orders, C, S, vals, grad, A44d, dA44d)
    elseif d == 3
        eval_UC_spline_3d(a, b, orders, C, S, vals, grad, A44d, dA44d)
    elseif d == 4
        eval_UC_spline_4d(a, b, orders, C, S, vals, grad, A44d, dA44d)
    end

    return (vals, grad)

end

# problem with this approach: the functions don't get cached.

for d = 1:4
    eval(create_function(d,"natural"))
end

for d = 1:4
    eval(create_function_with_gradient(d,"natural"))
end
