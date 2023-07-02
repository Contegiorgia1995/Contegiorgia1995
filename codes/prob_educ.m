function [Ve,Ce,Se,V0,tau] = prob_educ(par,options,Vp,Cp)
% Solve education choice

% Age 16 (age of HS) 
Ve        = cell(3,3);
Ce        = cell(3,3);
Se        = cell(3,3);

j_pos     = par.Je1_pos;
S         = par.grids{1,j_pos};
FE_pos     = par.inc.fe_pos;
PSY       = par.psy_val_hs;

V0        = zeros(length(S),length(FE_pos),length(PSY));
tau       = zeros(length(S),length(FE_pos),length(PSY));

%% Case 1: HS Dropout
educ = 1;
[Vhsd,Veaux,Ceaux,Seaux] = prob_educ_hs_drop(par,options,Vp,Cp); %#ok<*ASGLU>
for jp = 1:3
    Ve{educ,jp} = Veaux{1,jp};
    Ce{educ,jp} = Ceaux{1,jp};
    Se{educ,jp} = Seaux{1,jp};
end

%% Case 2: HS Graduate
educ = 2;
[Veaux,Ceaux,Seaux] = prob_educ_hs_grad(par,options,Vp,Cp);
for jp = 1:3
    Ve{educ,jp} = Veaux{1,jp};
    Ce{educ,jp} = Ceaux{1,jp};
    Se{educ,jp} = Seaux{1,jp};
end

%% Case 3: College Graduate
educ = 3;
[Veaux,Ceaux,Seaux] = prob_educ_co_grad(par,options,Vp,Cp);
for jp = 1:3
    Ve{educ,jp} = Veaux{1,jp};
    Ce{educ,jp} = Ceaux{1,jp};
    Se{educ,jp} = Seaux{1,jp};
end

%% Optimal education:
for ipsy = 1:length(PSY)
    for is = 1:length(S)
        for ife = 1:length(FE_pos)
            v_hsd = Vhsd(is,ife);
            v_hsg = Ve{2,1}(is,ife,ipsy);
            v_cg  = Ve{3,1}(is,ife,ipsy);
            
            v_aux = [v_hsd, v_hsg, v_cg];
            [V0(is,ife,ipsy),tau(is,ife,ipsy)] = max(v_aux);
        end
    end
end


end