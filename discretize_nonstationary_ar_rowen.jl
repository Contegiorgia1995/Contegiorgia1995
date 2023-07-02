using LinearAlgebra

function discretize_nonstationary_ar_rouwen(rho, sigma_z_l,sigma_z_t, n)
    q = 1/2 * (1+rho*sigma_z_l/sigma_z_t)
    nu = sqrt(n-1) * sigma_z_t
    P = [q 1-q;1-q q]
    for i in 2:n-1
    
        P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) P;zeros(1,i+1)] + (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P]
        
        P[2:i,:] = P[2:i,:]/2
    end
    
    zgrid = LinRange(-nu, nu, n)
    zgrid = reshape(zgrid, (3, 1))
    


    return P, zgrid