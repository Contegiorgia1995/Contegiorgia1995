function [m_model_sorted,m_title_data,m_mod_id] = moments_model(par,pol,mu,mu_cs,mu_ige,mu_chetty,...
                                                Inc_par_grid,mu_educ_ige,dist_age,IGE,options)
                                            
%% Compute moments from the model
% June 2015
switch options.timer_on
    case {'Y'}
        fprintf('\n Moments Model \n \n');
end

%% Distribute policy functions
% pol_mat             = fieldnames(pol);
% for i = 1 : length(pol_mat)
%     eval([cell2mat(pol_mat(i)) '= pol.(cell2mat(pol_mat(i)));']);
% end

[m_data,~,~] = moments_data;
Nm              = size(m_data,1);

% Empty array
count_m   = 1;
m_model   = NaN(Nm,1);
m_mod_id  = NaN(Nm,1);
m_title   = cell(Nm,1);

%% 1. Mean and Variance of Income by Education and Age
[ count_m,m_model,m_title,m_mod_id,mean_income_educ ] = mom_income_educ( par,count_m,m_model,m_title,m_mod_id,mu_cs,dist_age );

%% 2. Mean and Variance of Income by Age
[ count_m,m_model,m_title,m_mod_id ] = mom_income_age( par,count_m,m_model,m_title,m_mod_id,mu_cs );

%% 3. Educational attainement
[ count_m,m_model,m_title,m_mod_id,dist_educ ] = mom_dist_educ( par,count_m,m_model,m_title,m_mod_id,mu );

%% 4. Mean Income
[ count_m,m_model,m_title,m_mod_id, me_inc ] = mom_me_inc( par,count_m,m_model,m_title,m_mod_id,mean_income_educ,dist_educ );

%% 5: Fertility: mean, mean by educ and pr(fert = x)
[ count_m,m_model,m_title,m_mod_id ] = mom_fertility( par,count_m,m_model,m_title,m_mod_id,mu,options );

%% 6: Transfers to children
[ count_m,m_model,m_title,m_mod_id ] = mom_transfers( par,count_m,m_model,m_title,m_mod_id,mu,me_inc,IGE );

%% 7a1: Fertility elasticity wrt income: deciles
[ count_m,m_model,m_title,m_mod_id ] = mom_fert_elast_labinc( par,count_m,m_model,m_title,m_mod_id,mu,pol );

%% 7a2: Fertility elasticity wrt income: deciles by educ
[ count_m,m_model,m_title,m_mod_id ] = mom_fert_elast_labinc_byeduc( par,count_m,m_model,m_title,m_mod_id,mu,pol );

%% 9. Intergeneration Mobility: ranks of 100 + Log Log regression + Theil-L index 1 
[ count_m,m_model,m_title,m_mod_id ] = mom_rank_rank( par,count_m,m_model,m_title,m_mod_id,mu_chetty,Inc_par_grid );
[ count_m,m_model,m_title,m_mod_id ] = mom_rank_rank_ige( par,count_m,m_model,m_title,m_mod_id,IGE);

%% 10 Theil-L index
% Types: Parents educ and income
% Outcome: income at Jc
% Distribution: [mu_chetty_large,Inc_par_grid]
[ count_m,m_model,m_title,m_mod_id ] = mom_theil_ParInc( par,count_m,m_model,m_title,m_mod_id,mu_chetty,Inc_par_grid,options);


% keyboard
[ count_m,m_model,m_title,m_mod_id ] = mom_bottom_educ_inc( par,count_m,m_model,m_title,m_mod_id,mu_chetty,Inc_par_grid,options);

% Types: initial conditions
% Outcome: average lifetime earnings
% Distribution: [mu_LE,Grid_Avg_LE]
% [ count_m,m_model,m_title,m_mod_id ] = mom_theil_LE( par,count_m,m_model,m_title,m_mod_id,mu_LE,Grid_Avg_LE,options );

% Types: initial conditions
% Outcome: Income at Jc
% Distribution: [mu_IC_Inc,grid_IC_Inc]
% [ count_m,m_model,m_title,m_mod_id ] = mom_theil_IC_Inc( par,count_m,m_model,m_title,m_mod_id,mu_IC_Inc,grid_IC_Inc,options );

% Types: Parents educ and income
% Outcome: average lifetime earnings
% Distribution: [mu_Par_LE,Inc_par_grid,Grid_Avg_LE]
% [ count_m,m_model,m_title,m_mod_id ] = mom_theil_Par_LE( par,count_m,m_model,m_title,m_mod_id,mu_Par_LE,Inc_par_grid,Grid_Avg_LE,options);

%% Educ Parents - Educ Children
[ count_m,m_model,m_title,m_mod_id ] = mom_educ_ige( par,count_m,m_model,m_title,m_mod_id,mu_educ_ige,options );

%% 13. Income Inequality: Whole population   
[ count_m,m_model,m_title,m_mod_id ] = mom_ineq_all( par,count_m,m_model,m_title,m_mod_id,mu_cs );

%% 14. Income Inequality: By Age
% keyboard
[ count_m,m_model,m_title,m_mod_id ] = mom_ineq_byage( par,count_m,m_model,m_title,m_mod_id,mu_cs );

%% 15. Retirement Moments: Gov't benefits, Transfers and Savings. By Educ group.
[ count_m,m_model,m_title,m_mod_id ] = mom_ret_byeduc( par,count_m,m_model,m_title,m_mod_id,mu,options );

%% 16. Retirement Moments: Gov't benefits, Transfers and Savings. By income groups
% [ count_m,m_model,m_title,m_mod_id ] = mom_ret_byinc( par,count_m,m_model,m_title,m_mod_id,mu,options );

%% 17. Childcare Moments
[ count_m,m_model,m_title,m_mod_id ] = mom_childcare( par,count_m,m_model,m_title,m_mod_id,mu,pol );

%% 18. Education Variance Decomposition
[ count_m,m_model,m_title,m_mod_id ] = mom_educ_var_decomp( par,count_m,m_model,m_title,m_mod_id,mu_ige,options,pol );

%% 18b. Initial H and parents' characteristics
% [ count_m,m_model,m_title,m_mod_id ] = mom_parents_initialh( par,count_m,m_model,m_title,m_mod_id,mu_ige0,options,pol );

%% 19. Lifetime Earnings per individual
[ count_m,m_model,m_title,m_mod_id ] = mom_LE( par,count_m,m_model,m_title,m_mod_id,IGE,options );

%% 19b. Lifetime Earnings per individual: Age 24
[ count_m,m_model,m_title,m_mod_id ] = mom_LE_age( par,count_m,m_model,m_title,m_mod_id,IGE,par.Je2_pos+2,options );

%% 20: Fertility elasticity wrt total income
[ count_m,m_model,m_title,m_mod_id ] = mom_fert_elast_inctot( par,count_m,m_model,m_title,m_mod_id,mu,pol );

%% 21: Aggregates 
[count_m,m_model,m_title,m_mod_id,Agg_K] = Aggregates(par,mu_cs,options,count_m,m_model,m_title,m_mod_id);

%% 22. Aggregate Consumption and Aggregate Income
[ count_m,m_model,m_title,m_mod_id ] = mom_agg( par,count_m,m_model,m_title,m_mod_id,mu_cs,pol,Agg_K,options );

%% 6: Borrowing in College
[ count_m,m_model,m_title,m_mod_id ] = mom_borr_col( par,count_m,m_model,m_title,m_mod_id,mu,pol.Se{3,3} );

%% 6: Borrowing  Adults
[ count_m,m_model,m_title,m_mod_id ] = mom_borr( par,count_m,m_model,m_title,m_mod_id,mu_cs );

%% 23. Returns to education
[ count_m,m_model,m_title,m_mod_id ] = mom_educ_ret( par,count_m,m_model,m_title,m_mod_id,mu_cs);

%% Print
[m_data,m_title_data,~] = moments_data;
Nm              = size(m_data,1);
m_model_sorted  = nan(Nm,1);
[~,pos] = sort(m_mod_id);
for i = 1:count_m
    m_model_sorted(m_mod_id(i)) = m_model(i);
end

switch options.print_moms 
    case 'Y'
        % Print results
        fprintf('\t\t\t------------ MOMENTS ------------\n');
        fprintf('\t\t\t                                         Data    Model \n');
        for i = 1:count_m
            fprintf('%3.0f. \t %s \t %3.3f \t %3.3f \n',m_mod_id(pos(i)),char(m_title_data{m_mod_id(pos(i))}),m_data(m_mod_id(pos(i))),m_model(pos(i)));
        end

end
end