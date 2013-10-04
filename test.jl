using splines


print("\n1D interpolation\n\n")
orders = [100]
smin = [0.0]
smax = [1.0]
N_tries = 10

f(x) = sin(x*10)
# 1d test
a = linspace(0,1,orders[1])
b = f(a)

di = 1
M = length(a)

coefs = filter_coeffs(smin,smax,orders,b)

#s = linspace(-0.1,1.1,100000)
s = linspace(0.01,0.99,100000)

interp_values = eval_UBspline(smin,smax,orders,coefs,s)

tic()
for i = 1:N_tries
    interp_values = eval_UBspline(smin,smax,orders,coefs,s)
end

true_values = f(s)

print("Maximum interpolation error : ", norm(true_values-interp_values,Inf),"\n")
toc()

#using PyPlot
#plot(s,sol)
#plot(s,true_values)
#plot(a,b,"o")
#side effect of pyplot: 0.0 not recognized as a constant

# 2d test
print("\n2D interpolation\n\n")
smin = [0.0,0.0]
smax = [1.0,1.0]
orders = [10,20]
a1 = linspace(0,1,orders[1])
a2 = linspace(0,1,orders[2])
b = [sin(i+j) for i=a1, j=a2]
b = convert(Array{Float64},b)
#b = reshape(b, orders[1], orders[2])
#di = [1.0,1.0]
coefs = filter_coeffs(smin,smax,orders,b)

forders = [100,100]
N = prod(forders)
mgrid = [ [x y] for x=linspace(0,1,forders[1]), y=linspace(0,1,forders[2])]

vgrid = vcat(mgrid...)

interp_vals = eval_UBspline(smin,smax,orders,coefs,vgrid)
true_vals = [sin(x+y) for x=linspace(0,1,forders[1]), y=linspace(0,1,forders[2])]
true_vals = convert(Array{Float64},true_vals[:])

print("Maximum interpolation errori : ", norm(true_vals-interp_vals,Inf),"\n")

#tic()
#toc()

#coefs = filter_coeffs_2d(di, b)
#tic()
#coefs = filter_coeffs_2d(di, b)
#toc()

#using PyCall

#@pyimport test as ptt


#using Winston
#p = FramedPlot()
#add(p, Curve(a, b, "color","red"))
#add(p, Curve(a, coefs))
#Winston.display(p)



