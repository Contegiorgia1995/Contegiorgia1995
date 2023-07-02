function [mu,mu_cs,dist_age] = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options)
% Find distribution for a given cohort
switch options.timer_on
    case {'Y'}
        fprintf('\n Find ergodic distribution for cohort \n \n');
end

%% distribute policy functions
pol_mat = fieldnames(pol_dense);
for i = 1 : length(pol_mat)
    eval([cell2mat(pol_mat(i)) '= pol_dense.(cell2mat(pol_mat(i)));']);
end
clear pol_dense

dist_age  = zeros(par2.Jd_pos-1,1);
mu_cs     =  cell(par2.Jd_pos-1,1);

%% Age j=Jc + 3: Jr-1: 
for j_pos = par2.Jc_pos+par2.Je1_pos:par2.Jr_pos-1
    [mu{j_pos},Q_ergo{j_pos-1}] = trans_work(par2,mu{j_pos-1},S{j_pos-1},j_pos,Q_ergo{j_pos-1});
    switch options.timer_on
        case {'Y'}
            tot = sum(mu{j_pos}(:));
            fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
    end
end

%% age j = Jr: draw of transfers
j_pos = par2.Jr_pos;
mu{j_pos} = trans_workT(par2,mu{j_pos-1},S{j_pos-1},j_pos);
switch options.timer_on
    case {'Y'}
        tot = sum(mu{j_pos}(:));
        fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
end


%% Age j=Jr + 1: state are transfers
for j_pos = par2.Jr_pos+1:par2.Jd_pos
    mu{j_pos} = trans_ret(par2,mu{j_pos-1},S{j_pos-1},j_pos);
    switch options.timer_on
        case {'Y'}
            tot = sum(mu{j_pos}(:));
            fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
    end
end

pop_g   = pop_final^(1/par2.Jc);
sum_mu  = 0;
for j_pos = 1:par2.Jd_pos-1
    yy        = par2.age(j_pos);
    mu_cs{j_pos} = mu{j_pos}/(pop_g^yy);
    sum_mu       = sum_mu + sum(mu_cs{j_pos}(:));
end

for j_pos = 1:par2.Jd_pos-1
    mu_cs{j_pos} = mu_cs{j_pos}/sum_mu;
    dist_age(j_pos) = sum(mu_cs{j_pos}(:));
end

end