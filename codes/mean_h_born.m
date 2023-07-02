function [ out ] = mean_h_born( par,mu )
FE          = exp(par.inc.fe{1,1});
j_pos       = par.Je1_pos - 1;
inc_prob    = sum(sum(mu{j_pos}(:,:,:),1),3);
inc_prob    = inc_prob./sum(inc_prob(:));

out         = inc_prob(:)'*FE(:);
    

end