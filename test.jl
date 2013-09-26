include("splines_filter.jl")

# 1d test
a = linspace(0,1,10)
b = sin(a)

di = 1

M = length(a)

coefs = filter_coeffs_1d(di, b)
display(coefs)

# 2d test
a1 = linspace(0,1,100)
a2 = linspace(0,1,50)
b = [sin(i+j) for i=a1, j=a2]
b = convert(Array{Float64},b)
b = reshape(b, 100, 50)
di = [1.0,1.0]

tic()
toc()

coefs = filter_coeffs_2d(di, b)
tic()
coefs = filter_coeffs_2d(di, b)
toc()

#using PyCall

#@pyimport test as ptt


#using Winston
#p = FramedPlot()
#add(p, Curve(a, b, "color","red"))
#add(p, Curve(a, coefs))
#Winston.display(p)



