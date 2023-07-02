function [TR_PAR_VAL,TR_PAR_PDF,max_val] = dist_trans(par,tau0,options)
% Compute endogenous distribution of OAS transfers
% keyboard
M      = par.Lt;
% N_quad = 10;
%% Step 1: Child income distribution conditional on child's FE and parent transfers phi
% ch_earn_val (phi,h) : all possible values given (phi,FE)
% ch_earn_prob(phi,h) : corresponding probabilities

% 1.1 Compute human capital increments until age of transfers
Np    = length(par.prob{1,1}(par.Je1_pos,:));
max_N = length(par.Je1_pos :par.Jt_pos-1);

tau_max = max(tau0(:));
TreeP   = zeros(tau_max,Np^max_N);
TreeH   = zeros(tau_max,Np^max_N);

for educ=1:tau_max
    h0  = 1;
    p0  = 1;
    for j_pos = par.Je1_pos + (educ-1):par.Jt_pos-1
        p1 = par.prob{educ,1}(j_pos,:);
        h1 = par.deltas{educ,1}(j_pos,:);
        
        h0 = kron(h0,h1);
        p0 = kron(p0,p1);
    end
    TreeP(educ,1:length(p0)) = p0;
    TreeH(educ,1:length(h0)) = h0;
end

% 1.2 Obtain human capital after education
ch_earn_val  = zeros(length(par.PHI),length(par.H0),length(par.psy_val_hs),Np^max_N);
ch_earn_prob = zeros(length(par.PHI),length(par.H0),length(par.psy_val_hs),Np^max_N);
for ipsy = 1:length(par.psy_val_hs)
    for ih0=1:length(par.H0)
        for ip=1:length(par.PHI)
            educ = tau0(ip,ih0,ipsy);
            h   = par.H0(ih0);
            t   = 1;
            while t<educ
                h = par.Reducs(t,1) * h.^par.Reducs(t,2);
                t = t+1;
            end
            ch_earn_val(ip,ih0,ipsy,:)  = h.*TreeH(educ,:);
            ch_earn_prob(ip,ih0,ipsy,:) = TreeP(educ,:);
        end
    end
end

%% Step 2: Child income distribution conditional on parent transfers phi, parent educ and parents group
% ch_earn_val_2 (phi,EDUC,group,:) : all possible values given (phi,EDUC,FE)
% ch_earn_prob_2(phi,EDUC,group,:) : corresponding probabilities
EDUC           = par.educ;
GROUP          = par.cutoffs;

N_opt          = Np^max_N;
ch_earn_val_2  = zeros(length(par.PHI),length(EDUC),length(GROUP),N_opt*length(EDUC));
ch_earn_prob_2 = zeros(length(par.PHI),length(EDUC),length(GROUP),N_opt*length(EDUC));
aux_val        = nan(N_opt*length(EDUC),1); %This is educ of child
aux_prob       = nan(N_opt*length(EDUC),1); %This is educ of child
m_max          = 1;
for iphi = 1:length(par.PHI)
    for educ_p  = 1:length(EDUC)
        for ig_p = 1:length(GROUP)
            h0_prob = par.PG(educ_p,ig_p,:);
            pos = 1;
            for ipsy = 1:length(par.psy_val_hs)
                for ih0 = 1:length(par.H0)
                    aux_val(pos:pos+N_opt-1)                    = squeeze(ch_earn_val(iphi,ih0,ipsy,:));
                    aux_prob(pos:pos+N_opt-1)                   = h0_prob(ih0) * par.psy_prob(educ_p,ipsy) .* squeeze(ch_earn_prob(iphi,ih0,ipsy,:));
                    pos = pos+N_opt;
                    %                 ch_earn_val_2(iphi,educ_p,pos:pos+N_opt-1)  = squeeze(ch_earn_val(iphi,ih0,ipsy,:));
                    %                 ch_earn_prob_2(iphi,educ_p,pos:pos+N_opt-1) = h0_prob(ih0) * par.psy_prob(ipsy) .* squeeze(ch_earn_prob(iphi,ih0,ipsy,:));
                end
            end
            [val,~,ic] = unique(aux_val);
            m_aux       = length(val);
            ch_earn_val_2(iphi,educ_p,ig_p,1:m_aux)  = val;
            ch_earn_prob_2(iphi,educ_p,ig_p,1:m_aux) = accumarray(ic,aux_prob);
            m_max = max(m_max,m_aux);
        end
    end
end

switch options.Fertility
    case 'Endo'
        N_max = length(par.N);
        N_vec = par.N;
        n_1   = 2;
    case 'Exo'
        N_max = par.N;
        N_vec = 1:N_max;
        n_1   = 1;
end

PHI         = par.PHI;
TR_PAR_VAL  = cell(length(EDUC),N_max,length(PHI),length(GROUP));
TR_PAR_PDF  = cell(length(EDUC),N_max,length(PHI),length(GROUP));
curv        = 1.5;
max_val     = 0;
for iphi = 1:length(PHI)
    for educ = 1:length(EDUC)
        for iN = 1:N_max
            N = N_vec(iN);
            for ig = 1:length(GROUP)
                if N == 0
                    TR_PAR_VAL{educ,iN,iphi,ig}  = 0;
                    TR_PAR_PDF{educ,iN,iphi,ig}  = 1;
                elseif N == 1
                    aux_val                   = par.w .* par.xi .* (par.fam_size/2) .* squeeze(ch_earn_val_2(iphi,educ,ig,:));
                    aux_pdf                   = squeeze(ch_earn_prob_2(iphi,educ,ig,:));
                    [aux_val,pos]             = sort(aux_val);
                    aux_pdf                   = aux_pdf(pos);
                    TR_PAR_VAL{educ,iN,iphi,ig}  = ((linspace( min(aux_val(:))^(1/curv),(max(aux_val(:)))^(1/curv),M)).^curv)';
                    fspace_trans              = fundef({'spli', full(TR_PAR_VAL{educ,iN,iphi,ig}),0,1});
                    Q_tau                     = funbas(fspace_trans,aux_val);
                    TR_PAR_PDF{educ,iN,iphi,ig}  = (aux_pdf'* Q_tau)';
                    
                elseif N>1
                    % S(n) = sum(transfer of child i=1,...,n)
                    % S(n) = S(n-1) + tr(n)
                    % we know the distribution of S(n) and tr(n)!
                    
                    % Create grid of transfers for (iN,ip,ie,ife):
                    aux                          = N .* par.w .* par.xi .* (par.fam_size/2) .* squeeze(ch_earn_val_2(iphi,educ,ig,:));
                    TR_PAR_VAL{educ,iN,iphi,ig}  = ((linspace( min(aux(:))^(1/curv),(max(aux(:)))^(1/curv),M)).^curv)';
                    %                     TR_PAR_PDF{fe,educ,iN,iphi}  = zeros(1,M)';
                    
                    fspace_trans     = fundef({'spli', TR_PAR_VAL{educ,iN,iphi,ig},0,1});
                    
                    % S(n-1) distribution:
                    tr_val_prev_aux           = TR_PAR_VAL{educ,iN-1,iphi,ig};
                    tr_val_prev_aux           = squeeze(tr_val_prev_aux);
                    tr_pdf_prev_aux           = TR_PAR_PDF{educ,iN-1,iphi,ig};
                    tr_pdf_prev_aux           = squeeze(tr_pdf_prev_aux );
                    
                    % New draw:
                    tr_val_temp               = sparse(squeeze(TR_PAR_VAL{educ,n_1,iphi,ig}));
                    tr_pdf_temp               = sparse(squeeze(TR_PAR_PDF{educ,n_1,iphi,ig}));
                    
                    % Possible values of S(n)
                    E0                        = ones(1,length(tr_val_prev_aux));
                    E1                        = ones(1,length(tr_val_temp));
                    tr_val_aux                = kron(tr_val_temp',E0)+kron(E1,tr_val_prev_aux');
                    tr_prob_aux               = kron(tr_pdf_temp',E0).*kron(E1,tr_pdf_prev_aux');
                    
                    Q_tau                     = funbas(fspace_trans,tr_val_aux');
                    %                 TR_PAR_PDF{educ,iN,iphi}  = min((tr_prob_aux* Q_tau),0)'; %Min because sometimes get -1e-20...
                    TR_PAR_PDF{educ,iN,iphi,ig}  = max((tr_prob_aux* Q_tau),0)'; %max because sometimes get -1e-20...
                end % if iN
                max_val = max(max_val,max(TR_PAR_VAL{educ,iN,iphi,ig}(:)));
                %                 figure(1)
                %                 plot(TR_PAR_VAL{iN,ip,ie,ife},cumsum(TR_PAR_PDF{iN,ip,ie,ife}))
                %                 tit = strcat('N = ',num2str(iN),'phi = ',num2str(ip),'educ = ',num2str(ie),' FE = ',num2str(ife));
                %                 title(tit)
                %                 keyboard
            end
        end % for iN
    end % for ie
end % for ip




switch options.Fertility
    case 'Exo'
        TR_PAR_VAL  = TR_PAR_VAL(:,end,:,:);
        TR_PAR_PDF  = TR_PAR_PDF(:,end,:,:);
end


end