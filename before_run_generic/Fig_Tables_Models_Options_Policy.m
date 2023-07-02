% Figures and Tables from Exercise Model Options
clear; close all; clc; fignr = 1;
ex_version       	= 6;

%% Paths and results
addpath([pwd,'/run_generic']);
results_fold        = [pwd,'/results/'];

res_file_name       = strcat('results_model_options_policy_data_2000_',num2str(ex_version),'.mat');
file_model          = strcat(results_fold,'\',res_file_name);
res_2000            = load(file_model);

%% Data
[m_data_2000,mtitles,inc]       = moments_data;
[m_data_1960,~,~]               = moments_data_1960;
[m_data_fert_1960,~,~]          = moments_data_fert_1960;

% mtitles{284}    = 'Variance of lifetime earnings explained by initial conditions';

mtitles{108}       = '\hspace{.5in} Children income rank';
mtitles{109}       = '\hspace{.5in} Children dropouts';
mtitles{110}       = '\hspace{.5in} Children High school gradutes';
mtitles{111}       = '\hspace{.5in} Children college graduates';


moments_change_scale            = [3:5 109:111];
res_2000.res_MOM(:,moments_change_scale) = 100* res_2000.res_MOM(:,moments_change_scale);
m_data_2000(moments_change_scale) = 100* m_data_2000(moments_change_scale);
res_2000.res_MOM(:,98)      = res_2000.res_MOM(:,98).*res_2000.res_MOM(:,1);
res_ini_assets              = res_2000.res_MOM(:,98) + [ 0 0 20000]';
%% Tables fertility and mobility: levels & changes, long table
model_vec = 1;

for i_model = model_vec 
    if i_model == 1
        res = res_2000;
        ex_name = '2000_v2';
    end

    % Table 1: Levels  - bench - constant fert - constant transfers
    selected_models                 = [1 2 3];
    fid = fopen(['Table_policy_level_1b_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\multicolumn{1}{c}{\\textbf{Fertility Transfer}} & \\textbf{Initial Transfer} \\\\ \n');
    fprintf(fid,'\\hline\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Fertility and transfers}} \\\\ \n');
    selected_mom_1                  = [6 7 98 99];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n','Initial assets',res_ini_assets);
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Mobility}} \\\\ \n');
    
    fprintf(fid,' \\multicolumn{4}{l}{Bottom parents (income Q1 \\& dropouts)}\\\\ \n');
    
    selected_mom_1                  = [065 108:111];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    
    selected_mom_1                  = [3:5 14 121];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Inequality}} \\\\ \n');
    selected_mom_1                  = [124 284:287 289:290  56:60];
    mtitles{283}        = 'CV of Lifetime Earnings';
    mtitles{284}        = '\hspace{.5in} \% expl. by all initial conditions';
    mtitles{285}        = '\hspace{1in} \% expl. by human capital';
    mtitles{286}        = '\hspace{1in} \% expl. by transfers';
    mtitles{287}        = '\hspace{1in} \% expl. by school taste';
    mtitles{289}        = '\hspace{.5in} \% expl. by adult income shocks';
    mtitles{290}        = '\hspace{.5in} Total';
    res.res_MOM(:,289)    = 100 - res.res_MOM(:,284);
    res.res_MOM(:,290)    = 100;
    
    
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
    
    
        selected_models                 = [1 2 3];
    fid = fopen(['Table_policy_level_1_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\multicolumn{1}{c}{\\textbf{Fertility Transfer}} & \\textbf{Initial Transfer} \\\\ \n');
    fprintf(fid,'\\hline\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Fertility and transfers}} \\\\ \n');
    selected_mom_1                  = [6 7 98 99];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Mobility}} \\\\ \n');
    
    fprintf(fid,' \\multicolumn{4}{l}{Bottom parents (income Q1 \\& dropouts)}\\\\ \n');
    
    selected_mom_1                  = [065 108:111];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    
    selected_mom_1                  = [3:5 14 121];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Inequality}} \\\\ \n');
    selected_mom_1                  = [124 284:287 289:290  56:60];
    res.res_MOM(:,289)    = 100 - res.res_MOM(:,284);
    res.res_MOM(:,290)    = 100;
    
    
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
    
    % Table 2: changes  - bench - constant fert - constant transfers
    
    selected_models                 = [1 2 3 ];
    fid = fopen(['Table_policy_changes_1_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\multicolumn{1}{c}{\\textbf{Fertility Transfer}} & \\textbf{Initial Transfer} \\\\ \n');
    fprintf(fid,' & &  \\multicolumn{2}{c}{(change from benchmark)}  \\\\ \n');
    fprintf(fid,'\\hline\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Fertility and transfers}} \\\\ \n');
    selected_mom_1                  = [6 7 98 99];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
        fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Mobility}} \\\\ \n');
    
    fprintf(fid,' \\multicolumn{4}{l}{Bottom parents (income Q1 \\& dropouts)}\\\\ \n');
    
    selected_mom_1                  = [065 108:111];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
        fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
    end
    fprintf(fid,' && & \\\\ \n');
    
    selected_mom_1                  = [3:5 14 121];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
        fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');  
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Inequality}} \\\\ \n');


    selected_mom_1                  = [124 284:287 289:290  56:60];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        res_models = [res.res_MOM(selected_models(1),mom); (100*res.res_MOM(selected_models(2:end),mom)./res.res_MOM(selected_models(1),mom)-100)];
        fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
    end
    
    mom = 284;
    res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
    fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);


    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
end

%% Tables fertility and mobility: Slides
res = res_2000;
% ex_name = '2000';

selected_models                 = [1 2 3];
fid = fopen(['Table_policy_changes_1_',ex_name,'_slides.tex'],'w');
fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\multicolumn{1}{c}{\\textbf{Fertility Transfer}} & \\textbf{Initial Transfer} \\\\ \n');
fprintf(fid,' & &  \\multicolumn{2}{c}{(change from benchmark)}  \\\\ \n');
fprintf(fid,'\\hline\\hline\n');



fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Mobility}} \\\\ \n');

fprintf(fid,' \\multicolumn{4}{l}{Bottom parents (income Q1 \\& dropouts)}\\\\ \n');
mtitles{108}       = '\hspace{.5in} Children income rank';
mtitles{109}       = '\hspace{.5in} Children dropouts';
mtitles{110}       = '\hspace{.5in} Children High school gradutes';
mtitles{111}       = '\hspace{.5in} Children college graduates';

selected_mom_1                  = [108:111];
for n_mom = 1:length(selected_mom_1)
    mom = selected_mom_1(n_mom);
    res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
    fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
end
fprintf(fid,' && & \\\\ \n');

selected_mom_1                  = [3:5 14 121];
for n_mom = 1:length(selected_mom_1)
    mom = selected_mom_1(n_mom);
    res_models = [res.res_MOM(selected_models(1),mom); res.res_MOM(selected_models(2:end),mom) - res.res_MOM(selected_models(1),mom)];
    fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
end
fprintf(fid,' && & \\\\ \n');
fprintf(fid,'\\hline\n');


fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Inequality}} \\\\ \n');

mom = 124;
res_models = [res.res_MOM(selected_models(1),mom); (100*res.res_MOM(selected_models(2:end),mom)./res.res_MOM(selected_models(1),mom)-100)];
fprintf(fid,' %s   & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res_models);
fprintf(fid,' && & \\\\ \n');
fprintf(fid,'\\hline\n');


fprintf(fid,'\\end{tabular}');
fclose(fid);

