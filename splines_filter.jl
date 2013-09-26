function find_coefs_1d(delta_inv, M, data)

    # args: double, int, double[:], double[:

    N = M + 2
    
    data = data[:]

    if length(data) == N
        rhs = data
    else
        rhs = [0, data, 0]
    end

    basis = [1.0/6.0, 2.0/3.0, 1.0/6.0]

    vals = repeat( basis, outer=[M])
    co_x = repeat( [2:M+1], inner=[3])
    co_y = repeat( [1:3] , outer=[M]) + co_x - 2

    db = 4
    initial = [1.0,-2.0,1.0,0.0]*delta_inv*delta_inv
    final = [0.0,1.0,-2.0,1.0]*delta_inv*delta_inv

    vals = [initial, vals, final]
    co_x = [ones(Int,db), co_x, ones(Int,db)*(M+2)]
    co_y = [1:db, co_y, M+3-db:M+2]

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
        coefs[:,iy] = find_coefs_1d(dinv[1], Mx, data[:,iy])
    end

    # Now, solve in the Y-direction
    for ix in 1:Nx
        coefs[ix,:] =find_coefs_1d(dinv[2], My, coefs[ix,:])
    end

    return coefs

end

#function filter_coeffs_3d(double[:] dinv, double[:,:,:] data):

#    cdef int Mx = data.shape[0]
#    cdef int My = data.shape[1]
#    cdef int Mz = data.shape[2]

#    cdef int Nx = Mx+2
#    cdef int Ny = My+2
#    cdef int Nz = Mz+2

#    cdef double [:,:,:] coefs = np.zeros((Nx,Ny,Nz))

#    cdef int iy, ix, iz

#    # First, solve in the X-direction
#    for iy in range(My):
#        for iz in range(Mz):
#            find_coefs_1d(dinv[0], Mx, data[:,iy,iz], coefs[:,iy,iz])

#    # Now, solve in the Y-direction
#    for ix in range(Nx):
#        for iz in range(Mz):
#            find_coefs_1d(dinv[1], My, coefs[ix,:,iz], coefs[ix,:,iz])

#    # Now, solve in the Z-direction
#    for ix in range(Nx):
#        for iy in range(Ny):
#            find_coefs_1d(dinv[2], Mz, coefs[ix,iy,:], coefs[ix,iy,:])

#    return coefs

#end
