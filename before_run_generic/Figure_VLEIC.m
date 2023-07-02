% Figures and Tables from Exercise Model Options
clear; close all; clc; fignr = 1;
save_figs 			= 'Y';                        % If 'Y': save figures

%% Update the following paths:
addpath([pwd,'/run_generic']);
results_fold     = [pwd,'/results/'];

%% Path with codes
parent1 = fileparts(pwd);
parent1 = fileparts(parent1);
parent1 = strcat(parent1,'\',run_generic_path);
addpath(strcat(parent1,'\run_generic'));

%% Load results
file_model          = [fileparts(pwd),'\',results_fold,'\','LE_dist_age_111121'];
res                 = load(file_model);

%%
par                 = res.par;
start_age           = par.Je1_pos;
end_age             = par.Jr_pos-1;
steps               = 2;

age                 = par.age(start_age):steps:par.age(end_age);
age2                = par.age(start_age):8:par.age(end_age);

vle_ic_age          = res.VLE_IC_Age(start_age:end_age);

figure(fignr)
fig  = figure(fignr);
s1 = plot(age,vle_ic_age);
set(s1,'LineWidth',4)
xlabel('Age')
% ylabel('\fontsize{20} % of future earnings explained by current states')
set(gca,'XTick',age2)
box off;
axis([age(1) age(end) min(vle_ic_age) 100]);
set(gca,'FontSize',12)

graphname = 'fig_vle_ic_age';
print('-depsc2', sprintf('%s.eps',graphname));
saveas(fig, sprintf('%s.fig',graphname));
fignr  = fignr+1;