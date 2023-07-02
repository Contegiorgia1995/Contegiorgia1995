function rep = ret_rep(par,fe,educ,options)
% keyboard
% Pension replacement rate: 
% Source: Krueger, Ludwig 2015
% Parameters
w                   = par.w;
AGE_PROF            = par.inc.age_prof;
FE                  = par.inc.fe;

%% 1. Average lifetime income (AIME)
age_ret             = par.Je2_pos+1:par.Jr_pos-1;
totinc              = w * sum(exp(repmat(FE{educ,1}(fe)',length(age_ret),1) + repmat(AGE_PROF{educ,1}(age_ret),1,length(fe))));
    
prob_educ           = par.dist_educ;
den                 = 0;
h0_pos              = ceil(par.N_fe/2); % To find individual that starts with h0 = average
for ie = 1:length(par.educ)
    den             = den + prob_educ(ie) * w *sum(exp(FE{educ,1}(h0_pos)' + AGE_PROF{ie,1}(age_ret)));
end
y                   = totinc/den;

%% 2. Marginal replacement rates
switch options.ExPolRetirement
    case 'Y'
        tau1 = options.PolRetMult(1);
        tau2 = options.PolRetMult(2);
        tau3 = options.PolRetMult(3);
    case 'N'
        tau1   = 0.9;
        tau2   = 0.32;
        tau3   = 0.15;    
end


%% 3. Bendpoints
b1                  = 0.24;
b2                  = 1.35;
b3                  = 1.99;

%% 4. Replacement rate
rep_rate            = (tau1*y)                                .* (y<= b1) ... 
                    + (tau1*b1 + tau2*(y-b1))                 .* (y>b1) .* ( y<= b2) ...
                    + (tau1*b1 + tau2*(b2-b1) + tau3*(y- b2)) .* (y>b2) .* ( y<= b3) ...
                    + (tau1*b1 + tau2*(b2-b1) + tau3*(b3-b2)) .* (y>b3);
                
%% Replacement Benefits:
avg_inc             = den/length(age_ret);
rep                 = avg_inc .* rep_rate;
%         fprintf('educ = %i, fe = %i, y=%3.3f,ret_rep = %3.3f , rep = %3.3f, avg inc = %3.3f \n',educ,fe,y,rep_rate,rep,totinc/length(age_ret));
end