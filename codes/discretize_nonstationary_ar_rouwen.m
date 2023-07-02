function [zgrid, P] = discretize_nonstationary_ar_rouwen(rho, sigma_z_l,sigma_z_t, n)
% Rouwen method to discretize AR(1), non-stationary, process w symmetric innovations
% Based on WP Fella, Gallipoli and Pan (2015)
% It matches persistence of the original process exactly
    
    q = 1/2 * (1+rho*sigma_z_l/sigma_z_t);
    nu = sqrt(n-1) * sigma_z_t;

    P = [q 1-q;1-q q];

    for i = 2: n-1
        P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) P;zeros(1,i+1)] + ...
           (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P];
        P(2:i,:) = P(2:i,:)/2;
    end

    zgrid = linspace(-nu, nu, n)';

end