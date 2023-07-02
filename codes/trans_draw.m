function [mup, Q_ergo_out] = trans_draw(par,mu0,tau0,Q_ergo_in)
% From mu0 to mup
% keyboard
% Choice: education
j_pos     = par.Je1_pos;
S         = par.grids{1,j_pos};
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
PSY       = par.psy_val_hs;


%% 1. Choice of educ
mu0int          = zeros(length(S),length(FE_pos),length(PSY),length(EDUC));
mu0int(:,:,:,1) = mu0;    % initial distribution in extended grid, position 1
mu0int_vec      = mu0int(:);

switch Q_ergo_in.created 
    case 'N'
        % extend policy for education tau
        tau0p         = ones(length(S),length(FE_pos),length(PSY),length(EDUC));
        tau0p(:,:,:,1)  = tau0;
        tau0int_vec   = tau0p(:);

        % create transition matrix of educ: linear interpolant for policy of education
        fspaceergeduc   = fundef({'spli',EDUC,0,1});
        Qeduc           = funbas(fspaceergeduc,tau0int_vec);

        % create transition matrix of fixed states: Qfixed
        Qfixed              = kron(ones(length(EDUC),1),speye(length(S)*length(FE_pos)*length(PSY))); 

        % create aggregate transition: Q
        Q   = dprod(Qeduc,Qfixed); 
        Q_ergo_out.mat1 = Q;
        Q_ergo_out.created = 'Y';
    case 'Y'
        Q   = Q_ergo_in.mat1;
        Q_ergo_out = Q_ergo_in;
end
% new distribution: 
mu_int_vec      = Q' * mu0int_vec;
muint             = reshape(mu_int_vec,length(S),length(FE_pos),length(PSY),length(EDUC));

%% 2. draw of innovation z0
% load distribution of z0 conditional on education
z_prob     = par.inc.z_prob;
z0_prob    = [squeeze(z_prob{1,1}(par.Je1_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+2,1,:))'];      % (educ,prob)
INNO_pos   = par.inc.inno_pos;
NZ0        = length(INNO_pos);         

% vectorize distribution and policy
mu0int_vec    = muint(:);

switch Q_ergo_in.created 
    case 'N'
        % create transition matrix of Z0 conditional on education
        Qz                = kron(z0_prob,ones(length(S)*length(FE_pos)*length(PSY),1));

        % create transition matrix of fixed states: Qfixed
        Qfixed              = speye(length(S)*length(FE_pos)*length(PSY)*length(EDUC)); 

        % create aggregate transition: Q
        Q   = dprod(Qz,Qfixed); 
        Q_ergo_out.mat2 = Q;
        Q_ergo_out.created = 'Y';

    case 'Y'
        Q   = Q_ergo_in.mat2;
        Q_ergo_out = Q_ergo_in;
end

% new distribution: 
mup_vec = Q' * mu0int_vec;
mup     = reshape(mup_vec,length(S),length(FE_pos),length(PSY),length(EDUC),NZ0);

end