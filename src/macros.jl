function create_Phi(d, extrap)
    d = 1
    for i=1:d
        block = []
        rhs_1 = Symbol(string("tp_", i,"_",1))
        rhs_2 = Symbol(string("tp_", i,"_",2))
        rhs_4 = Symbol(string("tp_", i,"_",4))
        rhs_3 = Symbol(string("tp_", i,"_",3))
        for j=1:4
            lhs = Symbol(string("Phi_", i,"_",j))
            eq = :($lhs = (Ad[1,1]*rhs_1 + Ad[1,2]*rhs_2 + Ad[1,3]*rhs_3 + Ad[1,4]*rhs_4) )
            push!(block, eq)
        end
        println(block)
    end
end

create_Phi(3, "natural")


# def print_expr(symbs, inds=[]):
#     if len(symbs) == 0:
#         return 'C[{}]'.format(str.join(',',['i{}+{}'.format(i,k) for i,k in enumerate(inds)]))
#     else:
#         h = symbs[0]
#         q = symbs[1:]
#         exprs = [  '{}_{}*({})'.format(h,i,print_expr(q,inds + [i])) for i in range(4)]
#         return str.join( ' + ', exprs )

function tensor_prod(symbs, inds)
    if length(symbs)==0
        # return parse(string("C",inds,""))
        return (parse(string("C",string(inds),"")))
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
            println("e", e)
            push!(exprs,e)
        end
        return :(+$exprs)
    end
end

println(tensor_prod([], [1,2,3,4]))

println(tensor_prod(["a_1","a_2"], []))
/
