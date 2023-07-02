function [ c_final,sp_final,boundgrid ] = EGM_withchild( par,ucp, dispinc,Spgrid,S,n )
% EGM
r_sav  = par.r_sav;
r_debt = par.r_debt;
gammac = par.gammac;

% % Check ucp decreasing
% Lag =(ucp(1:end-1)<=1e+10) .* (ucp(2:end) - ucp(1:end-1));
% Lagmax = max(Lag);
% if Lagmax > 0 
%     fprintf('ucp is not decreasing\n')
% end

% Solve for a endogenous
c = ( ucp).^(-1/gammac);

lambdan = par.lambdan;
gamman  = par.gamman;
f_n     = (lambdan/n^(1-gamman))^(1/gammac);
c_k     = f_n * c;
a_endo  = (c + n*c_k + Spgrid' - dispinc );
a_endo = a_endo/(1+r_sav).*(a_endo>=0) + a_endo/(1+r_debt).*(a_endo<0);

% Check borrowing limits
a_bc    = a_endo(1); % Any current level of savings below this should be constrained
S_bc    = S(S < a_bc); % Borrowing constrained
c_bc    = ((1+r_sav) .* S_bc .* (S_bc >= 0) + (1+r_debt) .* S_bc .* (S_bc < 0) + dispinc - Spgrid(1))/(1+n*f_n);
c_k_bc  = f_n * c_bc;

% Interpolate in originial grid S
% (Since the borrowing limit is defined on strict inequality, all a_endo
% should be included in interpolation)

c_final    = approx_2d([S_bc'; a_endo], [c_bc'; c], S);
c_k_final  = approx_2d([S_bc'; a_endo], [c_k_bc'; c_k], S);
sp_final = (1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0) + dispinc - c_final - n*c_k_final; 
sp_final = sp_final .* ((sp_final>=Spgrid(1)) +(sp_final<Spgrid(1)-1e-6)) + ...
    Spgrid(1) .*(sp_final<Spgrid(1)).*(sp_final>=Spgrid(1)-1e-6);

if sum((sp_final<Spgrid(1)))>0
    fprintf('EGM error: savings are below debt limit - wrong extrapolation, min(s''): %3.2f \n',min(sp_final));
end

boundgrid = sum(sp_final>1.025*Spgrid(end));

% check for bad extrapolation
c_max   = ((1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0)+dispinc-Spgrid(1))/(1+n*f_n);
c_k_max   = f_n * c_max;
c_final = ( c_final<= c_max) .* c_final ...
    +     (c_final > c_max) .* c_max;
c_k_final = ( c_k_final<= c_k_max) .* c_k_final ...
    +     (c_k_final > c_k_max) .* c_k_max;

sp_final = (1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0) + dispinc - c_final - n*c_k_final; 

end