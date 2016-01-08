using Base.getindex
using Base.Cartesian



U(s,i) = Symbol(string(s,i))
U(s,i,j) = Symbol(string(s,i,"_",j))

function create_Phi(d, extrap)
    lines = []
    for i=1:d
        block = []
        rhs_1 = U("tp", i,1)
        rhs_2 = U("tp", i,2)
        rhs_4 = U("tp", i,4)
        rhs_3 = U("tp", i,3)
        if extrap == "none"
            for j=1:4
                eq = :($(U("Phi_",i,j)) = (Ad[$j,1]*$rhs_1 + Ad[$j,2]*$rhs_2 + Ad[$j,3]*$rhs_3 + Ad[$j,4]*$rhs_4) )
                push!(lines,eq)
            end
        elseif extrap == "natural"
            # block1 = []
            # for j=1:4
            #     eq = :($(U("Phi_",i,j)) = (Ad[1,1]*$rhs_1 + Ad[1,2]*$rhs_2 + Ad[1,3]*$rhs_3 + Ad[1,4]*$rhs_4) )
            #     push!(block1,eq)
            # end
            # block2 = []
            # for j=1:4
            #     eq = :($(U("Phi_",i,j)) = (Ad[1,1]*$rhs_1 + Ad[1,2]*$rhs_2 + Ad[1,3]*$rhs_3 + Ad[1,4]*$rhs_4) )
            #     push!(block2,eq)
            # end
            # block3 = []
            block = quote
                if $(U("t",i))<0
                    $( [ :($(U("Phi_",i,j)) = (dAd[$j,4]*$rhs_3 + Ad[$j,4]) ) for j=1:4 ]...)
                elseif $(U("t",i))>1
                    $( [ :($(U("Phi_",i,j)) = (3*Ad[$j,1] + 2*Ad[$j,2] + Ad[$j,3])*($rhs_3-1) + (Ad[$j,1]+Ad[$j,2]+Ad[$j,3]+Ad[$j,4]) ) for j=1:4 ]...)
                else
                    $( [ :($(U("Phi_",i,j)) = (Ad[$j,1]*$rhs_1 + Ad[$j,2]*$rhs_2 + Ad[$j,3]*$rhs_3 + Ad[$j,4]*$rhs_4) ) for j=1:4 ]...)

                end
            end
            # for l in block.args
            push!(lines,block.args...)

        end
    end
    return lines
end

create_Phi(3, "natural")

function tensor_prod(symbs, inds)
    if length(symbs)==0
        # return parse(string("C",inds,""))
        subscripts = [:($(U("i",i))+$(inds[i])) for i=1:length(inds)]
        return :(C[$(subscripts...)])
        # return (parse(string("C",string(inds),"")))
    else
        h = symbs[1]
        if length(symbs)>1
            q = symbs[2:end]
        else
            q = []
        end
        exprs = []
        for i = 1:4
            e = parse( string(h,"_",i,"*(",tensor_prod(q,cat(1,inds,[i])),")") )
            push!(exprs,e)
        end
        return :(+($(exprs...)))
    end
end

tensor_prod([], Int64[1,2,3,4])           # C[1,2,3,4]
tensor_prod(["Phi_1"], Int64[])  # Phi_1_1 * C[i1 + 1] + Phi_1_2 * C[i1 + 2] + Phi_1_3 * C[i1 + 3] + Phi_1_4 * C[i1 + 4]
tensor_prod(["Phi_1", "Phi_2"], Int64[]) # -> Phi_1_1 * (Phi_2_1 * C[i1 + 1,i2 + 1] + Phi_2_2 * C[i1 + 1,i2 + 2] + Phi_2_3 * C[i1 + 1,i2 + 3] + Phi_2_4 * C[i1 + 1,i2 + 4]) + Phi_1_2 * (Phi_2_1 * C[i1 + 2,i2 + 1] + Phi_2_2 * C[i1 + 2,i2 + 2] + Phi_2_3 * C[i1 + 2,i2 + 3] + Phi_2_4 * C[i1 + 2,i2 + 4]) + Phi_1_3 * (Phi_2_1 * C[i1 + 3,i2 + 1] + Phi_2_2 * C[i1 + 3,i2 + 2] + Phi_2_3 * C[i1 + 3,i2 + 3] + Phi_2_4 * C[i1 + 3,i2 + 4]) + Phi_1_4 * (Phi_2_1 * C[i1 + 4,i2 + 1] + Phi_2_2 * C[i1 + 4,i2 + 2] + Phi_2_3 * C[i1 + 4,i2 + 3] + Phi_2_4 * C[i1 + 4,i2 + 4])





function create_parameters(d)
    lines = []
    for i=1:d
        block = quote
            $(U("M",i)) = orders[$i]
            $(U("start",i)) =  a[$i]
            $(U("dinv",i)) = (orders[$i]-1.0)/(b[$i]-a[$i])
        end
        for ll in block.args
            push!(lines, ll)
        end
    end
    return lines
end

function create_local_parameters(d)
    lines = []
    for i=1:d
        bl = quote
            $(U("x",i)) = S[n,$i]
            $(U("u",i)) = ($(U("x",i)) - $(U("start",i)))*$(U("dinv",i))
            $(U("i",i)) = (floor(Int,$(U("u",i)) ))
            $(U("i",i)) = max( min($(U("i",i)),$(U("M",i))-2), 0 )
            $(U("t",i)) = $(U("u",i))-$(U("i",i))
            $(U("tp",i,1)) = $(U("t",i))*$(U("t",i))*$(U("t",i))
            $(U("tp",i,2)) = $(U("t",i))*$(U("t",i))
            $(U("tp",i,3)) = $(U("t",i))
            $(U("tp",i,4)) = 1.0;
        end
        for ll in bl.args
            push!(lines, ll)
        end
    end
    return lines
end

function create_function(d,extrap="natural")
    expr = quote
        function $(Symbol(string("eval_UC_spline_",d,"d")))( a, b, orders, C, S, V, Ad, dAd)
            $(create_parameters(d)...)
            N = size(S,1)
            # @simd( # doesn't seem to make any difference
            for n=1:N
                $(create_local_parameters(d)...)
                $(create_Phi(d,extrap)...)
                V[n] = $( tensor_prod([string("Phi_",i) for i=1:d], Int64[]) )
            end
            # )
        end
    end
    return expr
end

#
# @time fun = create_function(2,"natural");
# println(fun)
#
# d = 3
# K = 10
# N = 100000
# a = [0.0 for i in 1:d]
# b = [1.0 for i in 1:d]
# orders = [K for i in 1:d]
# C = rand([K+2 for i in 1:d]...)
# S = rand(N,d)*10-5
# V = zeros(N)
#
#
# Ad = [
#    [-1.0/6.0  3.0/6.0 -3.0/6.0 1.0/6.0];
#    [ 3.0/6.0 -6.0/6.0  0.0/6.0 4.0/6.0];
#    [-3.0/6.0  3.0/6.0  3.0/6.0 1.0/6.0];
#    [ 1.0/6.0  0.0/6.0  0.0/6.0 0.0/6.0]
# ]
#
# dAd = [
#    [ 0.0 -0.5  1.0 -0.5];
#    [ 0.0  1.5 -2.0  0.0];
#    [ 0.0 -1.5  1.0  0.5];
#    [ 0.0  0.5  0.0  0.0]
# ]
# #
# # d2A44d = [
# #    [ 0.0 0.0 -1.0  1.0];
# #    [ 0.0 0.0  3.0 -2.0];
# #    [ 0.0 0.0 -3.0  1.0];
# #    [ 0.0 0.0  1.0  0.0]
# # ]
#
# import splines
# using splines
# size(V)
#
# d = 2
# V0 = zeros(N)
# @time splines.eval_UC_spline_2d( a, b, orders, C, S, V0, Ad, dAd)
#
# V1 = zeros(N)
# vv1 = create_function(d, "natural")
# eval(vv1)
# eval_UC_spline_2( a, b, orders, C, S, V1, Ad, dAd)
# @time eval_UC_spline_2( a, b, orders, C, S, V1, Ad, dAd)
#
# vv1
# vv2 = create_function(d, "none")
# V2 = zeros(N)
# eval(vv2)
# eval_UC_spline_2( a, b, orders, C, S, V2, Ad, dAd)
# @time eval_UC_spline_2( a, b, orders, C, S, V2, Ad, dAd)
# # create_function(2, "none")
# V2 - V1
#
# V2 - V0
#
# V1 - V0
# maximum(S)
# minimum(S)
#
#
# vv1
