% Figures and Tables from Exercise Model Options
clear; close all; clc; fignr = 1;
ex_version       	= 6;

%% Paths and results
addpath([pwd,'/run_generic']);
results_fold     = [pwd,'/results/'];

res_file_name       = strcat('results_model_options_data_2000_',num2str(ex_version),'.mat');
file_model          = strcat(results_fold,'\',res_file_name);
res_2000            = load(file_model);

%% Data
[m_data_2000,mtitles,inc]       = moments_data;


% mtitles{284}    = 'Variance of lifetime earnings explained by initial conditions';

mtitles{108}       = '\hspace{.5in} Children income rank';
mtitles{109}       = '\hspace{.5in} Children dropouts';
mtitles{110}       = '\hspace{.5in} Children High school gradutes';
mtitles{111}       = '\hspace{.5in} Children college graduates';


moments_change_scale                        = [3:5 109:111];
res_2000.res_MOM(:,moments_change_scale)    = 100* res_2000.res_MOM(:,moments_change_scale);
m_data_2000(moments_change_scale)           = 100* m_data_2000(moments_change_scale);
res_2000.res_MOM(:,98)                      = res_2000.res_MOM(:,98).*res_2000.res_MOM(:,1);

%% Tables fertility and mobility: levels & changes, long table
model_vec = 1;
for i_model = model_vec % 2000, fert_1960, 1960
    
    if i_model == 1
        res = res_2000;
        ex_name = '2000';
    end    
    
    selected_models                 = [1 2 3];
    fid = fopen(['Table_model_level_1_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\textbf{Constant fertility} & \\textbf{Constant transfers} \\\\ \n');
    fprintf(fid,'\\hline\\hline\n');
    
    fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Fertility and transfers}} \\\\ \n');
    selected_mom_1                  = [6 7 98 99];
    for n_mom = 1:length(selected_mom_1)
        mom = selected_mom_1(n_mom);
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
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
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
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
        fprintf(fid,' %s & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
    end
    fprintf(fid,' && & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
    
    % Table 2: changes  - bench - constant fert - constant transfers
    
    selected_models                 = [1 3 9];
    fid = fopen(['Table_model_changes_1_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lccc}\n');
    fprintf(fid,' & \\textbf{Benchmark} & \\textbf{Constant fertility} & \\textbf{Constant transfers} \\\\ \n');
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
    
    selected_mom_1                  = [124 283];
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
ex_name = '2000';

selected_models                 = [1 3 9];
fid = fopen(['Table_model_changes_1_',ex_name,'_slides.tex'],'w');
fprintf(fid,'\\begin{tabular}{lccc}\n');
fprintf(fid,' & \\textbf{Benchmark} & \\textbf{Constant fertility} & \\textbf{Constant transfers} \\\\ \n');
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


%% Tables: H0, risk, school taste

res                     = res_2000;
ex_name                 = '2000';
selected_models         = [1 4:6 7 10:11];
fid = fopen(['Table_model_level_2_',ex_name,'.tex'],'w');
fprintf(fid,'\\begin{tabular}{lccccccc}\n');
fprintf(fid,' & \\textbf{Benchmark} & \\multicolumn{3}{c}{\\textbf{Initial human capital}} & \\textbf{No adult} & \\multicolumn{2}{c}{\\textbf{School taste}} \\\\ \n');
fprintf(fid,' & & \\textbf{constant} & \\textbf{iid} & \\textbf{no uncertainty}  & \\textbf{risk} &  \\textbf{iid} & \\textbf{constant} \\\\ \n');

fprintf(fid,'\\hline\\hline\n');

fprintf(fid,' \\multicolumn{8}{l}{\\textbf{Fertility and transfers}} \\\\ \n');
selected_mom_1                  = [6 7 98 99];
for n_mom = 1:length(selected_mom_1)
    mom = selected_mom_1(n_mom);
    fprintf(fid,' %s & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
end
fprintf(fid,' && & \\\\ \n');
fprintf(fid,'\\hline\n');

fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Mobility}} \\\\ \n');

fprintf(fid,' \\multicolumn{4}{l}{Bottom parents (income Q1 \\& dropouts)}\\\\ \n');

selected_mom_1                  = [65 108:111];
for n_mom = 1:length(selected_mom_1)
    mom = selected_mom_1(n_mom);
    fprintf(fid,' %s & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
end
fprintf(fid,' && & \\\\ \n');

selected_mom_1                  = [3:5 14 121];
for n_mom = 1:length(selected_mom_1)
    mom = selected_mom_1(n_mom);
    fprintf(fid,' %s & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f  \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
end
fprintf(fid,' && & \\\\ \n');
fprintf(fid,'\\hline\n');


fprintf(fid,' \\multicolumn{4}{l}{\\textbf{Inequality}} \\\\ \n');
selected_mom_1                  = [124 284:287 289:290 56:60];
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
    fprintf(fid,' %s & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f & %2.3f \\\\ \n',mtitles{mom},res.res_MOM(selected_models,mom));
end
fprintf(fid,' && & \\\\ \n');
fprintf(fid,'\\hline\n');



fprintf(fid,'\\end{tabular}');
fclose(fid);

%% Table decomposition

if do_only_2000 == 1
    model_vec = 1;
else
    model_vec = 1:2;
end

for i_model = model_vec % 2000, fert_1960, 1960
    
    if i_model == 1
        res = res_2000;
        ex_name = '2000';
    elseif i_model == 2
        res = res_fert_1960;
        ex_name = 'fert_1960';
    elseif i_model == 3
        res = res_1960;
        ex_name = '1960';
    end
    
    fid = fopen(['Table_model_change_decomp_',ex_name,'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lcc}\n');
    fprintf(fid,' & \\textbf{Mobility} & \\textbf{Inequality} \\\\ \n');
    fprintf(fid,' & {\\footnotesize Rank-rank correlation} & {\\footnotesize Var log(Life time earnings) } \\\\ \n');
    fprintf(fid,'\\hline\\hline\n');
    selected_mom_1                  = [121 124];

    fprintf(fid,' \\textbf{Benchmark} & %2.3f & %2.3f \\\\ \n',res.res_MOM(1,selected_mom_1));
    fprintf(fid,' & & \\\\ \n');
    fprintf(fid,'\\hline\n');
    
    
    fprintf(fid,' \\multicolumn{3}{l}{\\textbf{Exogenous}} \\\\ \n');
%     keyboard
    pos_model   = 4;
    res_models = [res.res_MOM(pos_model,selected_mom_1(1))- res.res_MOM(1,selected_mom_1(1))...
                 (100* res.res_MOM(pos_model,selected_mom_1(2))./res.res_MOM(1,selected_mom_1(2)) - 100)];
    fprintf(fid,' Constant initial human capital & %2.3f & %2.3f\\\\ \n',res_models);
    
    
    pos_model   = 7;
    res_models = [res.res_MOM(pos_model,selected_mom_1(1))- res.res_MOM(1,selected_mom_1(1))...
                 (100* res.res_MOM(pos_model,selected_mom_1(2))./res.res_MOM(1,selected_mom_1(2)) - 100)];
    fprintf(fid,' No adult income shocks & %2.3f & %2.3f\\\\ \n',res_models);
    
    pos_model   = 10;
    res_models = [res.res_MOM(pos_model,selected_mom_1(1))- res.res_MOM(1,selected_mom_1(1))...
                 (100* res.res_MOM(pos_model,selected_mom_1(2))./res.res_MOM(1,selected_mom_1(2)) - 100)];
    fprintf(fid,' Constant school taste & %2.3f & %2.3f\\\\ \n',res_models);
    
    
    fprintf(fid,' & & \\\\ \n');
    fprintf(fid,'\\hline\n');
    fprintf(fid,' \\multicolumn{3}{l}{\\textbf{Endogenous}} \\\\ \n');
    
    pos_model   = 3;
    res_models = [res.res_MOM(pos_model,selected_mom_1(1))- res.res_MOM(1,selected_mom_1(1))...
                 (100* res.res_MOM(pos_model,selected_mom_1(2))./res.res_MOM(1,selected_mom_1(2)) - 100)];
    fprintf(fid,' Constant fertility & %2.3f & %2.3f\\\\ \n',res_models);
    
    
    pos_model   = 8;
    res_models = [res.res_MOM(pos_model,selected_mom_1(1))- res.res_MOM(1,selected_mom_1(1))...
                 (100* res.res_MOM(pos_model,selected_mom_1(2))./res.res_MOM(1,selected_mom_1(2)) - 100)];
    fprintf(fid,' Constant transfers & %2.3f & %2.3f\\\\ \n',res_models);
    
    
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
    
end


%% Figure income inequality by age and models
age                     = 22:2:64;
% var_inc_m               = 100*(res_2000.res_MOM([3 8],883:902)-repmat(res_2000.res_MOM(1,883:902),2,1));
var_inc_m               = 100*(res_2000.res_MOM([3 9],883:904)./repmat(res_2000.res_MOM(1,883:904),2,1));


fig =figure(fignr);
f = plot(age,var_inc_m);
xlabel('Age')
ylabel('Income inequality (share of benchmark)')
set(f(1),'LineWidth',4,'LineStyle','-')
set(f(2),'LineWidth',4,'LineStyle','--')
ll=legend('Constant fertility','Constant transfers');
set(ll,'Location','NorthWest')
legend boxoff 
box off; grid off;
y_ticks         = [98 99 100];
set(gca,'YTick',y_ticks)
set(gca,'YTickLabel',sprintf('%1.0f %% \n',y_ticks))%             ax = gca;
axis([age(1) age(end) 97 100]) 
% axis tight
set(gca,'FontSize',12)

graphname = strcat('fig_inc_var_age_modeloptions');
print('-depsc2', sprintf('%s.eps',graphname));
saveas(fig, sprintf('%s.fig',graphname));
fignr  = fignr+1;

%% Transition matrix

if 1 == 1
    for n_model = [3 9]
        if n_model == 3
            fig_nam = 'fert_cte';
        elseif n_model == 9
            fig_nam = 'trans_cte';
        end
    
        pos_model           = 1;
        Q_model_bench      = zeros(5,5);
        Q_model_bench(1,:) = res_2000.res_MOM(pos_model,61:65);
        Q_model_bench(2,:) = res_2000.res_MOM(pos_model,66:70);
        Q_model_bench(3,:) = res_2000.res_MOM(pos_model,71:75);
        Q_model_bench(4,:) = res_2000.res_MOM(pos_model,76:80);
        Q_model_bench(5,:) = res_2000.res_MOM(pos_model,81:85);

        pos_model           = n_model;
        Q_model_fert_cte      = zeros(5,5);
        Q_model_fert_cte(1,:) = res_2000.res_MOM(pos_model,61:65);
        Q_model_fert_cte(2,:) = res_2000.res_MOM(pos_model,66:70);
        Q_model_fert_cte(3,:) = res_2000.res_MOM(pos_model,71:75);
        Q_model_fert_cte(4,:) = res_2000.res_MOM(pos_model,76:80);
        Q_model_fert_cte(5,:) = res_2000.res_MOM(pos_model,81:85);

        % Figure
        figure(fignr)
        fig = figure(fignr);    

        if 0 == 1
            hold on
            H_d = bar(Q_model_bench,'Grouped');
            set(H_d,'FaceColor','blue','LineStyle','none','BarWidth',.8);
            H_m = bar(Q_model_fert_cte,'Grouped');
            set(H_m,'FaceColor','red','LineStyle','none','BarWidth',.5);

            xlabel('\fontsize{30} Parent''s Quintile')
            hleglines = [H_d(1) H_m(1)];
            hleg = legend(hleglines,'Benchmark','Constant fertility');
            set(hleg,'Location','Best','Fontsize',30)
            legend boxoff
            axis tight
            box off; grid off;
            set(gca,'XTick',1:5)
        end
        
        if 1 == 1
            hold on
            H_d = bar(100*(Q_model_fert_cte(:,1)-Q_model_bench(:,1)),'Grouped');
            set(H_d,'FaceColor','blue','LineStyle','none','BarWidth',.8);
            
            xlabel('Parent''s Quintile')
            ylabel('\Delta Child''s Quintile Probability (p.p.) ')
            %         hleglines = [H_d(1) H_m(1)];
            %         hleg = legend(hleglines,'Benchmark','Constant fertility');
            %         set(hleg,'Location','Best','Fontsize',30)
            %         legend boxoff
            axis tight
            box off; grid off;
            %             set(gca,'YLim',[-3.52 3.0])
            set(gca,'YLim',[-4 4.0])
            set(gca,'XTick',1:5)
            set(gca,'FontSize',15)
            
            graphname = strcat('fig_transition_change_',fig_nam);
            print('-depsc2', sprintf('%s.eps',graphname));
            saveas(fig, sprintf('%s.fig',graphname));
            fignr  = fignr+1;
        end
        
        if 1 == 0
            hold on
            H_d = bar(100*(Q_model_fert_cte-Q_model_bench),'Grouped');
            set(H_d,'FaceColor','blue','LineStyle','none','BarWidth',.8);

            xlabel('Parent''s Quintile')
            ylabel('\Delta Child''s Quintile Probability (p.p.) ')
    %         hleglines = [H_d(1) H_m(1)];
    %         hleg = legend(hleglines,'Benchmark','Constant fertility');
    %         set(hleg,'Location','Best','Fontsize',30)
    %         legend boxoff
            axis tight
            box off; grid off;
%             set(gca,'YLim',[-3.52 3.0])
            set(gca,'YLim',[-4 4.0])
            set(gca,'XTick',1:5)
            set(gca,'FontSize',15)

        graphname = strcat('fig_transition_change_',fig_nam);
        print('-depsc2', sprintf('%s.eps',graphname));
        saveas(fig, sprintf('%s.fig',graphname));
        fignr  = fignr+1;
        end
        
        if 1 == 0
            hold on
            H_d = bar(100*(Q_model_bench),'Grouped');
            set(H_d,'FaceColor','blue','LineStyle','none','BarWidth',.8);
            
            xlabel('Parent''s Quintile')
            ylabel('Child''s Quintile Probability (p.p.) ')
            %         hleglines = [H_d(1) H_m(1)];
            %         hleg = legend(hleglines,'Benchmark','Constant fertility');
            %         set(hleg,'Location','Best','Fontsize',30)
            %         legend boxoff
            axis tight
            box off; grid off;
            set(gca,'YLim',[0 35])
            set(gca,'XTick',1:5)
            set(gca,'FontSize',15)
%             keyboard
            graphname = strcat('fig_transition_bench');
            print('-depsc2', sprintf('%s.eps',graphname));
            saveas(fig, sprintf('%s.fig',graphname));
            fignr  = fignr+1;
        end
    
        
    end
end