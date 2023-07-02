# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 18:27:11 2022

@author: Giorgia
"""

## Compute autocorrelation of income process in the model
import random

## 1. Simulate households
n_sim           = 10000
n_age           = (par.Jr_pos - 1)-(par.Je2_pos+2) + 1
sim_inc         = np.repeat(np.empty([3,n_sim]),n_age)
sim_inc         = np.reshape(sim_inc,[3,10000,22])
sim_inc[:]      = np.nan
shock_sim       = np.random.rand(3,n_sim,n_age)
shock_sim[:]    = random()

for educ in range (0,3):
    sim_inc[educ,:,1] = 1
    for age in range (2,n_age):        
        age_pos             = par.Je2_pos+age
        dist_cdf            = np.consum(par.prob{1}(age_pos,:));
        shock_val           = 1+sum((repmat(squeeze(shock_sim(educ,:,age)),length(dist_cdf),1) >= repmat(dist_cdf',1,n_sim)));
        sim_inc(educ,:,age) = sim_inc(educ,:,age-1).*par.deltas{educ,1}(age_pos,shock_val);
    end
end
sim_inc = log(sim_inc);

%% 2. Remove Age effects
% Create age
age_vec = repmat(par.age(par.Je2_pos+2:par.Jr_pos - 1),n_sim,1);
age2_vec = age_vec.^2;
% Matrix to save leftover
income_panel = nan(3,n_sim,n_age);
for educ =1:3
    incomes = squeeze(sim_inc(educ,:,:));
    Y       = incomes(:);
    X       = [ones(size(Y)) age_vec(:) age2_vec(:)];
    betas   = (X'*X)\(X'*Y);
    
    income_temp   = incomes - betas(2)*age_vec - betas(3)*age2_vec;
    income_panel(educ,:,:) = income_temp;
end

%% 3. Calculate Covariances
covs      = nan(3,n_age,n_age);
data_covs = nan(3,n_age*(n_age-1)/2,3); 
for educ = 1:3
    count = 1;
    for age = 1:n_age
        inc1 = income_panel(educ,:,age);
        for age2 = age:n_age
             inc2 = income_panel(educ,:,age2);
             aux                 = cov(inc1',inc2');
             covs(educ,age,age2) = aux(1,2);
             data_covs(educ,count,:) = [age age2-age covs(educ,age,age2)];
             count = count + 1;
        end
    end
end

%% 4. Calculate persistence
for educ = 1:3
    data = squeeze(data_covs(educ,:,:));
    
    x0 = [0.95 0.2 0.05 0.02];
    x0(1)   = log(x0(1)/(1 - x0(1)));
    x0(2:4) = log(x0(2:4));
    %x0(2:4)   = log(x0(2:4)./(0.3-x0(2:4)));
    f = @(x) moments_dif2(x,data);
    opt_fmin        = optimset('TolFun',1e-7,'MaxIter',10000,'MaxFunEvals',2000,'Display','off');
    [x1,~]          = fminunc(f,x0,opt_fmin);
    
    rho     = exp(x1(1))/(1+exp(x1(1)));
    var_z0  = exp(x1(2));
    var_m   = exp(x1(3));
    var_v   = exp(x1(4));
    
    fprintf('Education = %i \n rho = %g, var_z0 = %g, var_m = %g, var_v = %g \n',educ, rho, var_z0, var_m, var_v);
end