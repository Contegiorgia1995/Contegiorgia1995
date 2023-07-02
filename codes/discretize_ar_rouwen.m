function [zgrid, P] = discretize_ar_rouwen(mu_eps,rho, sigma_eps, n)
% Rouwen method to discretize AR(1) process w symmetric innovations
% It matches persistence of the original process exactly

%     mu_eps = 0;
    
    q = (rho+1)/2;
    nu = sqrt(n-1) * sigma_eps/(1-rho^2)^.5;

    P = [q 1-q;1-q q];

    for i = 2: n-1
        P = q*[P zeros(i,1);zeros(1,i+1)] + (1-q)*[zeros(i,1) P;zeros(1,i+1)] + ...
           (1-q)*[zeros(1,i+1); P zeros(i,1)] + q*[zeros(1,i+1); zeros(i,1) P];
        P(2:i,:) = P(2:i,:)/2;
    end

    zgrid = linspace(mu_eps/(1-rho)-nu, mu_eps/(1-rho)+nu, n)';

end