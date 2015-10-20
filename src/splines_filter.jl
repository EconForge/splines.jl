function filter_coeffs(smin::Vector{Float64}, smax::Vector{Float64}, orders::Vector{Int64}, data::Array{Float64})

    dinv = (smax - smin)./(orders.-1)
    d = length(smin)
    if d == 1
        coeffs = filter_coeffs_1d(dinv,data)
    elseif d == 2
        coeffs = filter_coeffs_2d(dinv,data)
    else
        error()
    end
    return coeffs
end

function find_coefs_1d(delta_inv, M, data)

    # args: double, int, double[:], double[:

    N = M + 2

    data = data[:]

    if length(data) == N
        rhs = data
    else
        rhs = [0; data; 0]
    end

    basis = [1.0/6.0, 2.0/3.0, 1.0/6.0]

    vals = repeat( basis, outer=[M])
    co_x = repeat( collect(2:M+1), inner=[3])
    co_y = repeat( collect(1:3) , outer=[M]) + co_x .- 2

    db = 4
    initial = [1.0,-2.0,1.0,0.0]*delta_inv*delta_inv
    final = [0.0,1.0,-2.0,1.0]*delta_inv*delta_inv

    vals = [initial; vals; final]
    co_x = [ones(Int,db); co_x; ones(Int,db)*(M+2)]
    co_y = [1:db; co_y; M+3-db:M+2]

    spmat = sparse(co_x, co_y, vals)

    sol = spmat \ rhs

    return sol


end

function filter_coeffs_1d(dinv, data)

  M = size(data,1)
  N = M+2

  coefs = find_coefs_1d(dinv[1], M, data)

  return coefs

end

function filter_coeffs_2d(dinv, data)

    Mx = size(data,1)
    My = size(data,2)

    Nx = Mx+2
    Ny = My+2

    coefs = zeros(Nx,Ny)


    # First, solve in the X-direction
    for iy in 1:My
        coefs[:,iy+1] = find_coefs_1d(dinv[1], Mx, data[:,iy])
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        coefs[ix,:] = find_coefs_1d(dinv[2], My, coefs[ix,:])
    end

    return coefs

end

function filter_coeffs_3d(dinv, data)

    Mx = size(data,1)
    My = size(data,2)
    Mz = size(data,3)

    Nx = Mx+2
    Ny = My+2
    Nz = Mz+2

    coefs = zeros(Nx,Ny,Nz)


    # First, solve in the X-direction
    for iy in 1:My
        for iz in 1:Mz
            coefs[:,iy+1,iz+1] = find_coefs_1d(dinv[1], Mx, data[:,iy,iz])
        end
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        for iz in 1:Mz
            ccoefs[ix,:,iz+1] = find_coefs_1d(dinv[2], My, coefs[ix,:,iz])
        end
    end

    # Now, solve in the Z-direction
    for ix in 1:Nx
        for iy in 1:Ny
            ccoefs[ix,iy,:] = find_coefs_1d(dinv[2], Mz, coefs[ix,iy,:])
        end
    end

    return coefs

end
