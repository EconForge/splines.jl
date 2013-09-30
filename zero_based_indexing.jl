
rewrite(x::Number) = x
rewrite(x::Symbol) = x
rewrite(x::String) = x
function rewrite(expr)
    if expr.head == :ref
        return Expr(expr.head, expr.args[1], [rewrite(:( 1 + $i)) for i in expr.args[2:end]]...)
    else
        return Expr(expr.head, [rewrite(i) for i in expr.args]...)
    end
end

macro zero_index(expr)
    eval( rewrite(expr) )
end

#@zero_index begin
#    a = [67,90]
#    b = [0,1]
#    print(a[0])
#    print(a[b[1]])
#end


