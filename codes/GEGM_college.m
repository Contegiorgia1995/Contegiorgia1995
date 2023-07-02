function [ c_final,sp_final,boundgrid ] = GEGM_college( par,ucp, dispinc,Spgrid,S,splVp,innop_prob )
% EGM
options = optimset('Display','off');

r_sav  = par.r_sav;
r_debt = par.r_debt;
gammac = par.gammac;
beta   = par.beta;
INNO_pos   = 1:length(innop_prob); %Because only in last period need to take expectations


%% Find no-concave region
difmin = [ucp(1:end-1) - ucp(2:end), 0];
indmin = (difmin<0);

difmax = [0, ucp(2:end) - ucp(1:end-1)];
indmax = (difmax>0);
if sum(indmin)>0
    vmin = min(ucp(indmin));
    vmax = max(ucp(indmax));
    %     vmin = min(ucp.*indmin+1e20.*(1-indmin));
    %     [vmax,~] = max(ucp.*indmax);
    
    [imin,~] = max((1:length(ucp)).*(ucp>vmax));
    
    if vmin > min(ucp)
        [~,imax] = min((1:length(ucp)).*(ucp<vmin)+1e20.*(ucp>=vmin));
    elseif vmin == min(ucp)
        imax = length(ucp);
        %     elseif vmin< min(ucp)
        %         keyboard
    end
    
else
    imin = 0;
    imax = 1;
end

% non concave region: [imin+1, imax-1]
% figure
% plot(Spgrid,ucp,Spgrid,vmin.*ones(size(Spgrid)),Spgrid,vmax.*ones(size(Spgrid)))



%% Solve for a endogenous
c = ( ucp').^(-1/gammac);
col_fact = (par.col_fact .* (Spgrid'<0) + 1 .* (Spgrid'>=0));
a_endo = (c + Spgrid'.*(1./col_fact) - dispinc );
a_endo  = a_endo/(1+r_sav).*(a_endo>=0) + a_endo/(1+r_debt).*(a_endo<0);

%
% figure(2)
% plot(a_endo,Spgrid)

%% Find global solution:
a_opt  = a_endo;
sp_opt = Spgrid;
c_opt  = c;

% Find global solution in non concave region
for i=imin+1:imax-1
    if i == imin + 1
        if length(innop_prob)>1
            [Sp3,inno3] = ndgrid(Spgrid,INNO_pos);
            Vpaux          = splVp(Sp3,inno3)*innop_prob';
        else
            Vpaux          = splVp(Spgrid)';
        end
        
        curv         = 2;
        Spgrid_dense_neg = -linspace(0,(-Spgrid(1))^(1/curv),300).^curv;
        Spgrid_dense = [Spgrid_dense_neg(end:-1:1) linspace(1e-6^(1/curv),max(Spgrid(end))^(1/curv),1000).^curv];
        col_fact_dense = (par.col_fact .* (Spgrid_dense'<0) + 1 .* (Spgrid_dense'>=0));
        
        
        if length(innop_prob)>1
            [Sp3,inno3] = ndgrid(Spgrid_dense,INNO_pos);
            Vpaux_dense = splVp(Sp3,inno3)*innop_prob';
        else
            Vpaux_dense = splVp(Spgrid_dense)';
        end
        
    end
    a_candidate = a_endo(i);
    
    c = (1+r_sav).*a_candidate.*(a_candidate>=0) + (1+r_debt).*a_candidate.*(a_candidate<0)  +dispinc - Spgrid'.*(1./col_fact);
    
    objfunc = c.^(1-gammac)/(1-gammac) .* (c>=0) ...
        + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(c<0) ...
        + beta*Vpaux;
    
    [~,amax] = max(objfunc);
    
    if amax ~= i % Discard solution
        %       fprintf('z:%3.2f,a''_candidate:%3.2f,a''_global:%3.2f, \n',a_candidate,Spgrid(i),Spgrid(amax))
        %         a_opt(i)  = nan;
        %         sp_opt(i) = nan;
        %         c_opt(i)  = nan;
        
        c = (1+r_sav).*a_candidate.*(a_candidate>=0) + (1+r_debt).*a_candidate.*(a_candidate<0)  +dispinc - Spgrid_dense'.*(1./col_fact_dense);
        
        objfunc = c.^(1-gammac)/(1-gammac) .* (c>=0) ...
            + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(c<0) ...
            + beta*Vpaux_dense;
        [~,amax] = max(objfunc);
        
        a_opt(i)  = a_candidate;
        sp_opt(i) = Spgrid_dense(amax);
        c_opt(i)  = c(amax);
    end
end

% figure(3)
% plot(a_endo,Spgrid,a_opt,sp_opt)

%% Borrowing constraint
% Case 1: Spgrid(1) is global solution
col_fact_1 = (par.col_fact .* (Spgrid(1)<0) + 1 .* (Spgrid(1)>=0));
if isnan(sp_opt(1)) == 0
    %     fprintf('Borrowing constraint - global \n')
    a_bc    = a_opt(1); % Any current level of savings below this should be constrained
    S_bc    = S(S < a_bc); % Borrowing constrained
    c_bc    = (1+r_sav) .* S_bc .* (S_bc >= 0) + (1+r_debt) .* S_bc .* (S_bc < 0) + dispinc - Spgrid(1)*(1/col_fact_1);
else % Approximate a_bc
    [~,pos] = min(sp_opt);
    a_opt1    = a_opt(pos);
    if a_opt1<S(1);
        %         fprintf('Borrowing constraint - local, but never binding \n')
        a_bc    = a_opt(pos);
        S_bc    = S(S < a_bc); % Borrowing constrained
        c_bc    = (1+r_sav) .* S_bc .* (S_bc >= 0) + (1+r_debt) .* S_bc .* (S_bc < 0) + dispinc - Spgrid(1)*(1/col_fact_1);
    else
        %         fprintf('Borrowing constraint - local & potentially binding \n')
        S_bc = S(S<a_opt1);
        
        
        if length(innop_prob)>1
            [Sp3,inno3] = ndgrid(Spgrid,INNO_pos);
            Vpaux          = splVp(Sp3,inno3)*innop_prob';
        else
            Vpaux          = splVp(Spgrid)';
        end
        
        c = repmat((1+r_sav) .* S_bc .* (S_bc >= 0) + (1+r_debt) .* S_bc .* (S_bc < 0),size(Spgrid')) +dispinc - repmat(Spgrid'.*(1./col_fact),size(S_bc));
        
        objfunc = c.^(1-gammac)/(1-gammac) .* (c>=0) ...
            + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(c<0) ...
            + beta*repmat(Vpaux,size(S_bc));
        
        [~,amax] = max(objfunc,[],1);
        Sp_vfi   = Spgrid(amax);
        col_fact_vfi = (par.col_fact .* (Sp_vfi<0) + 1 .* (Sp_vfi>=0));
        c_bc    = (1+r_sav) .* S_bc .* (S_bc >= 0) + (1+r_debt) .* S_bc .* (S_bc < 0) + dispinc - Sp_vfi.*(1./col_fact_vfi);
    end
end


% Discard local solutions and sort grid
a_opt = a_opt((isnan(a_opt)==0));
c_opt = c_opt((isnan(c_opt)==0));
[a_opt,ord] = sort(a_opt);
c_opt = c_opt(ord);

% Interpolate in originial grid S
% (Since the borrowing limit is defined on strict inequality, all a_endo
% should be included in interpolation)

c_final  = approx_2d([S_bc'; a_opt], [c_bc'; c_opt], S);
sp_final = (1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0) + dispinc - c_final;
col_fact = (par.col_fact .* (sp_final<0) + 1 .* (sp_final>=0));
sp_final = sp_final.*col_fact; % Back to original grid
sp_final = sp_final .* ((sp_final>=Spgrid(1)) +(sp_final<Spgrid(1)-1e-6)) + ...
    Spgrid(1) .*(sp_final<Spgrid(1)).*(sp_final>=Spgrid(1)-1e-6);


if length(innop_prob)>1
    [Sp3,inno3] = ndgrid(sp_final,INNO_pos);
    Vpaux0          = splVp(Sp3,inno3)*innop_prob';
else
    Vpaux0          = splVp(sp_final);
end

objfunc0 = c_final.^(1-gammac)/(1-gammac) .* (c_final>=0) ...
    + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(c_final<0) ...
    + beta*Vpaux0;

boundgrid = sum(sp_final>1.025*Spgrid(end));

% Check if there are slow decreases in C, i.e., not jumps
difc = [0; c_final(2:end) - c_final(1:end-1)];
inddecc = (difc<0);
inddecc = (inddecc.*(([1; inddecc(1:end-1)]==1) + ([inddecc(2:end); 1] == 1))>0);
inddecc(1:end-2) = max([inddecc(1:end-2)  inddecc(3:end)],[],2);
inddecc(3:end) = max([inddecc(3:end)  inddecc(1:end-2)],[],2);
inddecc(objfunc0<0) = 1;
inddecc([0; objfunc0(1:end-1)]<0) = 1;

grid_prob = 1:length(c_final);
do_interp = 1;
for i = grid_prob(inddecc)
    if do_interp == 1
        curv         = 2;
        Spgrid_dense_neg = -linspace(0,(-Spgrid(1))^(1/curv),300).^curv;
        Spgrid_dense = [Spgrid_dense_neg(end:-1:1) linspace(1e-6^(1/curv),max(Spgrid(end))^(1/curv),1000).^curv];
        col_fact_dense = (par.col_fact .* (Spgrid_dense'<0) + 1 .* (Spgrid_dense'>=0));
        
        
        if length(innop_prob)>1
            [Sp3,inno3] = ndgrid(Spgrid_dense,INNO_pos);
            Vpaux_dense = splVp(Sp3,inno3)*innop_prob';
        else
            Vpaux_dense = splVp(Spgrid_dense)';
        end
        
        do_interp    = 0;
    end
    
    
    c = (1+r_sav).*S(i).*(S(i)>=0) + (1+r_debt).*S(i).*(S(i)<0)  +dispinc - Spgrid_dense'.*(1./col_fact_dense);
    
    objfunc = c.^(1-gammac)/(1-gammac) .* (c>=0) ...
        + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(c<0) ...
        + beta*Vpaux_dense;
    [~,amax] = max(objfunc);
    
    sp_final(i) = Spgrid_dense(amax);
    c_final(i)  = c(amax);
end

% check for bad extrapolation
c_max   = (1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0)+dispinc-Spgrid(1)*(1/col_fact_1);
c_final = ( c_final<= c_max) .* c_final ...
    +     (c_final > c_max) .* c_max;

sp_final = (1+r_sav).*S'.*(S'>=0) + (1+r_debt).*S'.*(S'<0) + dispinc - c_final;
col_fact = (par.col_fact .* (sp_final<0) + 1 .* (sp_final>=0));
sp_final = sp_final.*col_fact; % Back to original grid

if ((min(sp_final)-Spgrid(1))/abs(Spgrid(1)) > 0.01) && (min(sp_final)-Spgrid(1)<0)
    fprintf('EGM College error: savings are below debt limit - wrong extrapolation, min(s''): %3.4f \n',100*min(sp_final-Spgrid(1))/abs(Spgrid(1)));
end
end